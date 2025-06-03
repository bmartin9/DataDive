#' @title sim_gaussian simulates a normal y from given data matrix X
#' @param X an n by p matrix
#' @param pve a scalar percentage variance explained
#' @param effect_num a scalar number of true nonzero effects
#' @return train_n a scalar number of trainning samples
#' @return sim_y an n vector simulated gaussian y
#' @return beta a p vector of effects
#' @return mean_corX mean of correlations of X (lower triangular entries of correlation matrix of X)
sim_gaussian = function(X, pve, effect_num, cv_idx = NULL, cv_cor = NULL){
    mafs = calc_maf(X)
    if(pve == 0.005){
        n = min(50000,dim(X)[1])
    }else if(pve == 0.02){
        n = min(12500,dim(X)[1])
    }else if( pve == 0.1){
        n = min(2500,dim(X)[1])
    }else if(pve == 0.3){
        n = min(800,dim(X)[1])
    }else{
        n = dim(X)[1]
    }

    X = center_scale(X[1:n,])
    p = dim(X)[2]

    if (is.null(cv_idx)){
        beta.idx = sample(p, effect_num)
    }else{
        beta.idx = cv_idx
    }

    beta = rep(0,p)
    beta.values = rnorm(effect_num)
    beta[beta.idx] = beta.values


    if (effect_num==1){
        mean_corX = 1
    } else {
        if (effect_num==0){
            mean_corX=NULL
            resid_var = 0
        }else{
            effectX = X[,beta.idx]
            corX = cor(effectX, use= "pairwise.complete.obs")
            mean_corX = mean(abs(corX[lower.tri(corX)]), na.rm = TRUE)
        }}

    if(effect_num==0){
        sigma = 1
        sim.y = rnorm(n, 0, 1)
        y = (sim.y - mean(sim.y))/sd(sim.y)
    } else {

        # take care of missing values
        X_na = X
        X_na[is.na(X)] = 0
        y_genetic = X_na %*% beta

        pheno_var = var(y_genetic) / pve
        resid_var = pheno_var - var(y_genetic)

        epsilon = rnorm(n, mean = 0, sd = sqrt(resid_var))
        y.sim = y_genetic + epsilon

        y = y.sim
        # beta = beta/sd(y.sim)
        # y = y.sim/sd(y.sim)
        # resid_var = resid_var/var(y.sim)
    }

    # calculate SumStats from values
    snp_ebt = snp_ese = NULL
    for (snp in 1:p){
        gwas_df = data.frame(Y = y, X = X[,snp]); # y is already centered, X is centered and scaled
        gwas_res = speedglm::speedglm(formula = Y~X,
                                      data = gwas_df,
                                      family = gaussian(link = 'identity'))
        summary_gwas_res = summary(gwas_res)
        snp_ebt = c(snp_ebt,summary_gwas_res$coefficients$Estimate[2])
        snp_ese = c(snp_ese,summary_gwas_res$coefficients$`Std. Error`[2])
    }

    return(list(Y = y, sigma2 = resid_var,
                beta = beta, mean_corX = mean_corX,
                betas = snp_ebt, ses = snp_ese,
                mafs = mafs, cv_idx = beta.idx, cv_n = effect_num,
                cv_cor = cv_cor, pve = pve, n=n))
}

# simulation functions
calc_maf = function(gg) {
    maf <- (colSums(gg == 1, na.rm=TRUE) + 2 * colSums(gg == 2, na.rm=TRUE)) / (2 * dim(gg)[1])
    return(maf)
}

center_scale = function(X){
    X = scale(X, center = TRUE, scale = TRUE)
    return(X)
}
