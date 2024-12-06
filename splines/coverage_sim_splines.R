# Load necessary libraries
library(mgcv)
library(rhdf5)
library(R.matlab)
library(foreach)
library(doParallel)

# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less core than available
cl <- makeCluster(num_cores)
registerDoParallel(cl)

bootstrap_moves_plus <- function(X, y, n_replicas, n_movements) {
  
  flag <- TRUE
  selected_i <- NULL
  
  # Repeat until selected indices have non-zero standard deviation
  while (flag) {
    selected_i <- sample(1:length(y), n_replicas, replace = TRUE)
    if (sd(y[selected_i]) != 0) {
      flag <- FALSE
    }
  }
  
  Xboot <- matrix(0, nrow = n_replicas, ncol = ncol(X))
  yboot <- numeric(n_replicas)
  
  # Iterate over selected indices
  for (i in 1:length(selected_i)) {
    index <- selected_i[i]
    x <- X[index, ]
    yi <- y[index]
    
    # Create x_movement
    x_movement <- unlist(lapply(seq_along(x), function(j) rep(j, x[j])))
    
    # Bootstrap x_movement
    x_movement_boot <- sample(x_movement, n_movements, replace = TRUE)
    
    # Initialize xboot
    xboot <- numeric(ncol(X))
    
    # Populate xboot with the sampled movements
    for (xk in x_movement_boot) {
      xboot[xk] <- xboot[xk] + 1
    }
    
    # Append to Xboot and yboot
    Xboot[i, ] <- xboot
    yboot[i] <- yi
  }
  
  return(list(Xboot = Xboot, yboot = yboot))
}


# Assume you have:
# - y: a vector of binary responses of length m
# - X: an m x n^2 design matrix
# - n: the dimension of your image (so n^2 pixels)
niter = 50
nboot = 48
nx = 30
n = 500
tn = 2880
n_new = 5000
tn_new = 720

for( iter in 1:niter){
  
  for(type in c("lake", "river", "sides")){
    
    filename = paste('sim_',type,'_bird_moves_n',n,'_snr1_nx',nx,'_ny',nx,'_tn',tn,'iter',iter,'.mat', sep = "")
    print(filename)
    
    boot_filename = paste('BOOT_sim_',type,'_bird_moves_n',n,'_snr1_nx',nx,'_ny',nx,'_tn',tn,'iter',iter,'nboot',nboot,'tn_new', tn_new, 'n_new', n_new,'.mat',sep="")
    spline_filename = paste('spline_',boot_filename, sep="")
    
    
    if( file.exists(spline_filename)){
      spline_data <- readMat(spline_filename)
      #if ("all.propspline.boot" %in% names(spline_data)) {
      #  next
      #}
    }
    
    X <- h5read(filename, "/X")
    y <- h5read(filename, "/y")
    
    #results = c()
    results <- foreach(iboot = 1:nboot, .packages = c("mgcv")) %dopar% {
    
    #for(iboot in 1:nboot){  
      set.seed(iboot)
      print(iboot)
      out=bootstrap_moves_plus(X, y, n_new, tn_new);
      

      # Create a data frame for modeling
      data <- data.frame(y = out$y)
      
      # Create image pixel coordinates
      pixel_grid <- expand.grid(x = 1:nx, z = 1:nx)
      
      # Since X is large, we can't include it directly in the model formula
      # Instead, we'll define a custom smooth over the image pixels
      
      # Define the penalty matrix for the image smoothness
      # Using a thin plate spline penalty for 2D smooth
      
      # Prepare the smooth term
      smooth_term <- smoothCon(s(x, z, bs = "tp", k = min(c(n_new,nx^2/10))), data = pixel_grid)[[1]]
      
      # Extract the penalty matrix
      S <- smooth_term$S[[1]]
      
      # Compute the smooth basis for each pixel
      X_smooth <- PredictMat(smooth_term, pixel_grid)
      
      # Now, compute the effective design matrix for the model
      # Multiply X (m x n^2) by X_smooth (n^2 x s), where s is the number of smooth basis functions
      # This gives us an m x s matrix
      
      X_effective <- out$Xboot %*% X_smooth
      
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
      
      fit <- gam(out$y~-1+X_effective, data = data, family = binomial(),
                 paraPen = paraPen_list)
      
      # Summary of the model
      summary(fit)
      
      # Extract the estimated image b
      # Compute b = X_smooth %*% coefficients
      b_estimated <- X_smooth %*% coef(fit)
      
      # Reshape b_estimated back to image form
      estimated_image <- matrix(b_estimated, nx, nx)
      #results[[iboot]] = estimated_image
      
    }
    
    all_propspline_boot = c()
    
    for(iboot in 1:nboot){
      all_propspline_boot = cbind(all_propspline_boot, results[[iboot]][1:nx^2])  
    }
    
    writeMat(spline_filename, all_propspline_boot=all_propspline_boot)
  }
}
stopCluster(cl)

