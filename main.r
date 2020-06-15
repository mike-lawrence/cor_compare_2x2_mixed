library(tidyverse)
library(rstan)
library(ezStan) #install via: remotes::install_github('mike-lawrence/ezStan') ; remotes::install_github('mike-lawrence/loggr')

#generate some fake data ----

#simulation parameters
N = 1000 #subject per group
cor_a1 = .1
cor_a2 = .6
cor_b1 = .4
cor_b2 = .9


#generate data
tibble(
	group = c('a','a','b','b')
	, condition = c('x','y','x','y')
	, cor = c(cor_a1,cor_a2,cor_b1,cor_b2)
) %>%
	dplyr::group_by(
		group
		, condition
	) %>%
	dplyr::summarise(
		as_tibble(
			MASS::mvrnorm(
				n = N
				, mu = c(0,0)
				, Sigma = matrix(c(1,cor,cor,1),2,2)
			)
		)
		, .groups = 'keep'
	) %>%
	#add a subject id (appending group bc each unique person should have a unique id)
	dplyr::mutate(
		id = factor(paste0(group[1],'_',1:n()))
	) ->
	obs_data

#check that the correlations are approximately as expected
obs_data %>%
	dplyr::group_by(
		group
		, condition
	) %>%
	dplyr::summarise(
		value = cor(V1,V2)
	)

# compute some transforms of the data ----

#reshape the data a bit
obs_data %>%
	dplyr::ungroup() %>%
	#gather and unite measure with condition for a 4-level variable
	tidyr::gather(
		key = measure
		, value = value
		, V1
		, V2
	) %>%
	tidyr::unite(
		measure_condition
		, measure
		, condition
	) %>%
	#spread to wide
	tidyr::spread(
		key = measure_condition
		, value = value
	) ->
	obs_data_wide

print(obs_data_wide)

#prep data for Stan
data_for_stan = list(
	nSubj = nrow(obs_data_wide)
	, nW = ncol(obs_data_wide)-2
	, nB = length(unique(obs_data_wide$group))
	, Y = obs_data_wide %>% dplyr::select(-group,-id)
	, groupMembership = as.numeric(factor(obs_data_wide$group)) #1/2 indicating group
)

# sample using Stan ----

#compile model
cor_compare_mixed = ezStan::build_stan('cor_compare_mixed.stan')

#start the sampling (automatically parallel by default)
ezStan::start_stan(
	data = data_for_stan
	, mod = cor_compare_mixed
)

#watch the sampling progress
ezStan::watch_stan()

#when done, collect
post = ezStan::collect_stan() #ignore warnings

check_hmc_diagnostics(post) #pay attention to any warnings here!

#check rhat and ess
monitor(post)

# #take a look at the summary for a few variables
# commented-out bc there's a bug in naming I just discovered
# ezStan::stan_summary(
# 	from_stan = post
# 	, par = 'coef_means'
# )
#
# ezStan::stan_summary(
# 	from_stan = post
# 	, par = 'coef_sds'
# )
#
# ezStan::stan_summary(
# 	from_stan = post
# 	, par = 'cor_mat'
# ) %>%
# 	print(n=Inf)

#get the posterior on cor_mat
cor_mat = rstan::extract(post,par='cor_mat')[[1]]
str(cor_mat) #first dim is sample number, 2nd is group, then 4x4 correlation matrix

group1_cor_diff = rep(NA,dim(cor_mat)[1])
group2_cor_diff = rep(NA,dim(cor_mat)[1])
var = obs_data_wide %>% dplyr::select(-group,-id) %>% names()
for(i in 1:dim(cor_mat)[1]){
	group1_cor_diff[i] =
		cor_mat[i,1,var=='V1_x',var=='V2_x'] -
		cor_mat[i,1,var=='V1_y',var=='V2_y']
	group2_cor_diff[i] =
		cor_mat[i,2,var=='V1_x',var=='V2_x'] -
		cor_mat[i,2,var=='V1_y',var=='V2_y']
}

cor_diffs = tibble(
	group1 = group1_cor_diff
	, group2 = group2_cor_diff
	, sample_num = 1:length(group1_cor_diff)
)

#view the posterior for each group's difference
cor_diffs %>%
	tidyr::gather(
		key = group
		, value = value
		, -sample_num
	) %>%
	ggplot()+
	facet_wrap(
		~ group
	)+
	geom_histogram(
		mapping = aes(
			x = value
		)
	)

#compute the difference-of-differences
cor_diffs %>%
	dplyr::mutate(
		double_diff = group1 - group2
	) ->
	cor_diffs

#histogram
cor_diffs %>%
	ggplot()+
	geom_histogram(
		mapping = aes(
			x = double_diff
		)
	)

#median & 95% credible interval
quantile(cor_diffs$double_diff,c(.5,.025,.975))
