library(tidyverse)
library(rstan)
library(ezStan) #install via: remotes::install_github('mike-lawrence/ezStan') ; remotes::install_github('mike-lawrence/loggr')

#generate some fake data ----

#simulation parameters
N = 500 #subject per group
cor_a1 = .5
cor_a2 = .5
cor_b1 = .2
cor_b2 = .6


#generate data
tibble(
	group = c('a','a','b','b')
	, condition = c('1','2','1','2')
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
	#ensure the data are sorted (I *think* the stan code assumes sorted-by-id-first)
	dplyr::arrange(
		id
		, condition
		, group
	) %>%
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
	) ->
	obs_data_long

#get the within-subject contrast matrix
obs_data_long %>%
	ezStan::get_contrast_matrix(
		formula = ~0+measure_condition #"0+" makes toggles indicator contrasts
		, contrast_kind = 'contr.treatment'
	) ->
	W

#double-check the unique contrasts are what we expect
W %>%
	tibble::as_tibble() %>%
	dplyr::distinct() %>%
	View()

#reduce data to 1 row per subject then generate between-Ss contrasts
obs_data_long %>%
	dplyr::group_by(
		id
		, group
	) %>%
	dplyr::summarize(
		.groups='drop'
	) ->
	obs_id_group

obs_id_group %>%
	get_contrast_matrix(
		formula = ~ 0 + group #"0+" makes toggles indicator contrasts
		, contrast_kind = 'contr.treatment'
	) ->
	B

#double-check the unique contrasts are what we expect
B %>%
	tibble::as_tibble() %>%
	dplyr::distinct() %>%
	View()

#prep data for Stan
data_for_stan = list(
	nSubj = length(unique(obs_data_long$id))
	, nY = nrow(W)
	, nW = ncol(W)
	, nB = ncol(B)
	, Y = obs_data_long$value
	, W = W
	, B = B
	, subjIndices = ezStan::get_subject_indices(obs_data_long$id)
	, groupMembership = as.numeric(factor(obs_id_group$group)) #1/2 indicating group
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
