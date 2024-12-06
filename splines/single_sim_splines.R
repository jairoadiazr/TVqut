# Load necessary libraries
library(mgcv)
library(rhdf5)
# Assume you have:
# - y: a vector of binary responses of length m
# - X: an m x n^2 design matrix
# - n: the dimension of your image (so n^2 pixels)
niter = 100

for( iter in 59:niter){
  
  for(nx in c(30,50)){
    
    for(tn in c(30, 30*24*4)){
      
      for(type in c("lake", "river", "sides")){
        
        for(n in c(500, 5000)){
          
          filename = paste('sim_',type,'_bird_moves_n',n,'_snr1_nx',nx,'_ny',nx,'_tn',tn,'iter',iter,'.mat', sep = "")
          print(filename)
          X <- h5read(filename, "/X")
          y <- h5read(filename, "/y")

          # Check if 'q' exists in the file
          if ("/Propspline" %in% h5ls(filename)$name) {
            # Delete variable 'q' if it exists
            h5delete(filename, "/Propspline")
            cat("Variable 'Propspline' was deleted.\n")
          } else {
            cat("Variable 'Propspline' does not exist in the file.\n")
          }
          
          # Create a data frame for modeling
          data <- data.frame(y = y)
          
          # Create image pixel coordinates
          pixel_grid <- expand.grid(x = 1:nx, z = 1:nx)
          
          # Since X is large, we can't include it directly in the model formula
          # Instead, we'll define a custom smooth over the image pixels
          
          # Define the penalty matrix for the image smoothness
          # Using a thin plate spline penalty for 2D smooth
          
          # Prepare the smooth term
          smooth_term <- smoothCon(s(x, z, bs = "tp", k = min(c(n,nx^2/5))), data = pixel_grid)[[1]]
          
          # Extract the penalty matrix
          S <- smooth_term$S[[1]]
          
          # Compute the smooth basis for each pixel
          X_smooth <- PredictMat(smooth_term, pixel_grid)
          
          # Now, compute the effective design matrix for the model
          # Multiply X (m x n^2) by X_smooth (n^2 x s), where s is the number of smooth basis functions
          # This gives us an m x s matrix
          
          X_effective <- X %*% X_smooth
          
          # Now, we can set up the GAM with the effective design matrix
          # We'll treat the smooth coefficients as the parameters to estimate
          
          # Create a data frame with the effective predictors
          for (j in 1:ncol(X_effective)) {
            data[[paste0("s", j)]] <- X_effective[, j]
          }
          
          # Create the formula
          smooth_vars <- paste0("s", 1:ncol(X_effective))
          formula <- as.formula(paste("y ~ -1 +", paste(smooth_vars, collapse = " + ")))
          
          # Set up the penalty for the smooth coefficients
          
          paraPen_list <- list(X_effective=list(rank=qr(S)$rank,S))
          data = c()
          data$X_effective = X_effective
          data$y = y
          # Fit the GAM
          
          fit <- gam(y~-1+X_effective, data = data, family = binomial(),
                     paraPen = paraPen_list)
          
          # Summary of the model
          summary(fit)
          
          # Extract the estimated image b
          # Compute b = X_smooth %*% coefficients
          b_estimated <- X_smooth %*% coef(fit)
          
          # Reshape b_estimated back to image form
          estimated_image <- matrix(b_estimated, nx, nx)
          h5write(estimated_image, filename, "/propspline")
        }
      }
    }
  }
}
