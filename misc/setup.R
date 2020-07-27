load("misc/grid.RData")

mu_0=0.001875
gamma=c(-2,runif(3))
longLat=c(535000,295000)
sigma=10000
kappa=-2.5

#st_make_grid does rows and columns weirdly
rows=30
cols=40
grid=st_make_grid(grid,n=c(cols,rows))
nTime=20

#st_make_grid indexing:
indices=t(matrix(((rows)*(cols)):1,ncol=rows))

#reorder grid:
grid=grid[c(indices)[1:((rows)*cols)]]


coords=grid%>%st_centroid()%>%st_coordinates()
X_reaction=cbind(matrix(1,nrow=nTime*(rows*cols),ncol=1),
                 matrix(runif(nTime*(rows*cols)*3),ncol=3))
diffusionType=0
lengthX=1
lengthY=1

save(mu_0,
     gamma,
     longLat,
     sigma,
     kappa,
     rows,
     cols,
     grid,
     nTime,
     coords,
     X_reaction,
     diffusionType,
     lengthX,
     lengthY,
     file="tests/testthat/setup.RData")
