test_that("diffusion works", {
  source("../../misc/diffusion.R")
  load("../../misc/setup.RData")
  
  for (diffusionType in -1:1){
    u0=computeDiffusion(mu_0,
                        alpha,
                        gamma,
                        longLat,
                        sigma,
                        kappa,
                        coords,
                        X_diffusion,
                        X_reaction,
                        rows,
                        cols,
                        nTime,
                        diffusionType,
                        TRUE,
                        lengthX,
                        lengthY,
                        TRUE)
    expect_equal(u0,u0)
  }
})
