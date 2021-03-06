
test_that("du_dgamma works", {
  source("../../misc/diffusion.R")
  load("../../misc/setup.RData")
  for (diffusionType in -1:1){
    du_dgamma0=du_dgamma(mu_0,
                         alpha,
                         gamma,
                         longLat,
                         sigma,
                         kappa,
                         coords,
                         X_diffusion,
                         X_reaction,
                         rows-2,
                         cols-2,
                         nTime,
                         diffusionType,
                         TRUE,
                         lengthX,
                         lengthY,
                         TRUE)
    
    step=0.00001
    for (i in 1:length(du_dgamma0)){
      gamma0=gamma
      gamma1=gamma
      gamma0[i]=gamma[i]-step
      gamma1[i]=gamma[i]+step
      u0=computeDiffusion(mu_0,
                          alpha,
                          gamma0,
                          longLat,
                          sigma,
                          kappa,
                          coords,
                          X_diffusion,
                          X_reaction,
                          rows-2,
                          cols-2,
                          nTime,
                          diffusionType,
                          TRUE,
                          lengthX,
                          lengthY,
                          TRUE)
      u1=computeDiffusion(mu_0,
                          alpha,
                          gamma1,
                          longLat,
                          sigma,
                          kappa,
                          coords,
                          X_diffusion,
                          X_reaction,
                          rows-2,
                          cols-2,
                          nTime,
                          diffusionType,
                          TRUE,
                          lengthX,
                          lengthY,
                          TRUE)
      du_dgamma1=(u1-u0)/(2*step)
      expect_equal(du_dgamma0[[i]],du_dgamma1,tol=5*sqrt(.Machine$double.eps))
    }
  }
})
