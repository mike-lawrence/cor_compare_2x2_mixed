functions{
	//a function to turn the long data to wide coefficients
	matrix get_coefs(matrix W, vector Y, int [,] Z) {
		matrix[size(Z),cols(W)] coefs ;
		for(i in 1:size(Z)){
			coefs[i] = transpose(
				inverse(
					transpose( W[Z[i,1]:Z[i,2]] )
					* W[Z[i,1]:Z[i,2]]
				)
				* transpose( W[Z[i,1]:Z[i,2]] )
				* Y[Z[i,1]:Z[i,2]]
			);
		}
		return coefs ;
	}
}
data {
	// nSubj: number of subjects
	int nSubj ;
	// nY: total number of observations
	int nY ;
	// nW: number of within-subject predictors
	int nW ;
	// nB: number of between-subject predictors
	int nB ;
	// Y: vector of observations
	vector[nY] Y ;
	// W: within-subject predictor matrix
	matrix[nY,nW] W ;
	// B: between-subject predictor matrix
	matrix[nSubj,nB] B ;
	// subjIndices: list of start and end indices in Y
	//  corresponding to data from each subject
	int subjIndices[nSubj,2] ;
	// groupMembership: list of each subject's group membership (hack)
	int groupMembership[nSubj] ;
}
transformed data{
	// scaled_coefs: matrix of wide-format coefficients
	//  for each subject given their data and predictors
	matrix[nSubj,nW] scaled_coefs ;
	//use the get_coefs() function, simultaneously scaling Y
	// for easy default priors
	scaled_coefs = get_coefs(W,(Y-mean(Y))/sd(Y),subjIndices) ;
}
parameters {
	// scaled_coef_means: z-scale coefficient means
	matrix[nB,nW] scaled_coef_means ;
	// scaled_coef_sds: z-scale coefficient sds (unique per level of B)
	vector<lower=0>[nW] scaled_coef_sds[nB] ;
	// cor_mat: correlations amongst coefficients (unique per level of B)
	corr_matrix[nW] cor_mat[nB] ;
}
model {
	// weakly informed priors
	to_vector(scaled_coef_means) ~ normal(0,1) ;
	for(b in 1:nB){
		scaled_coef_sds[b] ~ weibull(2,1) ; //peaked around .8; declines to zero at zero, not too heavy right tail
		cor_mat[b] ~ lkj_corr(2) ; //relatively (but not completely) flat prior on correlations
	}
	//coefficients as multivariate normal
	for(s in 1:nSubj){
		scaled_coefs[s,] ~ multi_normal(
			  B[s,] * scaled_coef_means
			, quad_form_diag(cor_mat[groupMembership[s]], scaled_coef_sds[groupMembership[s]])
		) ;
	}
}
generated quantities{
	// coef_means: coefficient means on the original
	//  scale of Y
	matrix[nB,nW] coef_means ;
	// coef_means: coefficient sds on the original
	//  scale of Y
	vector[nW] coef_sds[nB] ;
	for(b in 1:nB){
		coef_sds[b] = scaled_coef_sds[b] * sd(Y) ;
	}
	coef_means = scaled_coef_means * sd(Y) ;
	//adding back the mean to the intercept only
	coef_means[1,1] = coef_means[1,1] + mean(Y) ;
}
