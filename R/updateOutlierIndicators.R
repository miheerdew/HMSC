#' @importFrom stats dnorm pnorm dpois runif
updateOutlierIndicators = function(Y.orig,Beta,iSigma,Eta,Lambda, Loff,X,Pi,dfPi,distr,rL,odEps){

   ny = nrow(Y.orig)
   ns = ncol(Y.orig)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))

   switch(class(X)[1L],
          matrix = {
             LFix = X%*%Beta
          },
          list = {
             LFix = matrix(NA,ny,ns)
             for(j in 1:ns)
                LFix[,j] = X[[j]]%*%Beta[,j]
          }
   )
   LRan = vector("list", nr)
   for(r in seq_len(nr)){
      if(rL[[r]]$xDim == 0){
         LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
      } else{
         LRan[[r]] = matrix(0,ny,ns)
         for(k in 1:rL[[r]]$xDim)
            LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),k]) %*% Lambda[[r]][,,k]
      }
   }

   E = Reduce("+", c(list(LFix), LRan))
   if(!is.null(Loff)) E = E + Loff

   # Variable to store the log-likelihoods of each entry of Y
   likelihood_Y = matrix(NA,ny,ns)
   indNA = is.na(Y.orig)
   std = matrix(iSigma^-0.5,ny,ns,byrow=TRUE)

   # For species that follow a normal model
   indColNormal = (distr[,1]==1)
   likelihood_Y[,indColNormal] = dnorm(Y.orig[,indColNormal], mean=E[,indColNormal], sd=std[,indColNormal])

   # For species that follow a Probit model
   indColProbit = (distr[,1]==2)
   YProbit = as.logical(Y.orig[,indColProbit])
   EProbit = E[,indColProbit]
   # Use the Probit likelihood $\Phi(E)^{y}(1-\Phi(E))^{1-y}$
   likelihood_Y[,indColProbit] = ifelse(YProbit, pnorm(EProbit), 1 - pnorm(EProbit))

   # For species that follow a Poisson model
   indColPoisson = (distr[,1]==3)
   likelihood_Y[ , indColPoisson] = dpois(
            Y.orig[ , indColPoisson],
            lambda=exp(E[ , indColPoisson]))

   outlierProb = odEps/(odEps + (1-odEps)*likelihood_Y)
   outlierIndicators = matrix(runif(ny*ns) < outlierProb, ny, ns)

   return(outlierIndicators)
}
