data {
	// nSubj: number of subjects
	int nSubj ;
	// nW: number of within-subject predictors
	int nW ;
	// nB: number of between-subject predictors
	int nB ;
	// Y: matrix of observations in wide format (one row per subject)
	matrix[nSubj,nW] Y ;
	// groupMembership: list of each subject's group membership (hack)
	int groupMembership[nSubj] ;
}
transformed data{
	// scaled_Y: matrix of wide-format data
	matrix[nSubj,nW] scaled_Y ;
	// scaling Y for easy default priors
	scaled_Y = (Y-mean(Y))/sd(Y) ;
}
parameters {
	// scaled_Y_means: z-scale Y column means
	matrix[nB,nW] scaled_Y_means ;
	// scaled_coef_sds: z-scale Y column sds (unique per level of B)
	vector<lower=0>[nW] scaled_Y_sds[nB] ;
	// cor_mat: correlations amongst coefficients (unique per level of B)
	corr_matrix[nW] cor_mat[nB] ;
}
model {
	// weakly informed priors
	to_vector(scaled_Y_means) ~ normal(0,1) ;
	for(b in 1:nB){
		scaled_Y_sds[b] ~ weibull(2,1) ; //peaked around .8; declines to zero at zero, not too heavy right tail
		cor_mat[b] ~ lkj_corr(2) ; //relatively (but not completely) flat prior on correlations
	}
	//coefficients as multivariate normal
	for(s in 1:nSubj){
		scaled_Y[s,] ~ multi_normal(
			  scaled_Y_means[groupMembership[s]]
			, quad_form_diag(cor_mat[groupMembership[s]], scaled_Y_sds[groupMembership[s]])
		) ;
	}
}
generated quantities{
	// Y_means: column means on the original scale of Y
	matrix[nB,nW] Y_means ;
	// Y_means: coefficient sds on the original scale of Y
	vector[nW] Y_sds[nB] ;
	for(b in 1:nB){
		Y_sds[b] = scaled_Y_sds[b] * sd(Y) ;
	}
	Y_means = scaled_Y_means * sd(Y) + mean(Y) ;
}
