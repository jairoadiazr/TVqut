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

nboot = 100
weight = 1000
year = "ALL"
ov = ''
interpolation = 15

for(ehdv in 6){
  spline_filename = paste('spline_EHDV',ehdv,'.mat', sep="")
  print(spline_filename)
  
  points <- read.csv(paste0("points_EHDV-", ehdv, "_weight_", weight, "_year_", year, "_interpolation_", interpolation, "_", ov, ".csv"))
  points = as.matrix(points[,2:3])
  # Adjusting the points as per the calculation in MATLAB
  points[,1] <- round(points[,1] * weight - min(points[,1] * weight) + 1)
  points[,2] <- round(points[,2] * weight - min(points[,2] * weight) + 1)
  
  y <- read.csv(paste0("Y_vector_EHDV-", ehdv, "_weight_", weight, "_year_", year, "_interpolation_", interpolation, "_", ov, ".csv"))
  X <- read.csv(paste0("X_matrix_EHDV-", ehdv, "_weight_", weight, "_year_", year, "_interpolation_", interpolation, "_", ov, ".csv"))
  y = as.vector(y[,2])
  X = as.matrix(X[,-1])
  
  X_matrix = matrix(0,nrow=max(points[,1]),ncol=max(points[,2]))
  for(i in 1:nrow(points)){
    X_matrix[points[i,1], points[i,2]] = 1
  }
  
  results <- foreach(iboot = 1:nboot, .packages = c("mgcv")) %dopar% {
  #for(iboot in 1:nboot){
    
    n_replicas = 600
    n_movements = 2000
    
    out = bootstrap_moves_plus(X, y, n_replicas, n_movements)
    Xboot = out$Xboot
    yboot = out$yboot
    

    
    newX = c()
    for(j in 1:nrow(Xboot)){
      cur_X_vector = Xboot[j,]
      cur_X_matrix = matrix(0,nrow=max(points[,1]),ncol=max(points[,2]))
      for(i in 1:length(cur_X_vector)){
        cur_X_matrix[points[i,1], points[i,2]] = cur_X_vector[i]
      }
      newX = rbind(newX, cur_X_matrix[1:length(cur_X_matrix)])
    }
    
    # Create a data frame for modeling
    data <- data.frame(y = yboot)
    
    # Create image pixel coordinates
    nx = nrow(cur_X_matrix)
    ny = ncol(cur_X_matrix)
    
    pixel_grid <- expand.grid(x = 1:nx, z = 1:ny)
    
    # Since X is large, we can't include it directly in the model formula
    # Instead, we'll define a custom smooth over the image pixels
    
    # Define the penalty matrix for the image smoothness
    # Using a thin plate spline penalty for 2D smooth
    
    # Prepare the smooth term
    smooth_term <- smoothCon(s(x, z, bs = "tp", k = nx*ny/2), data = pixel_grid)[[1]]
    
    # Extract the penalty matrix
    S <- smooth_term$S[[1]]
    
    # Compute the smooth basis for each pixel
    X_smooth <- PredictMat(smooth_term, pixel_grid)
    
    # Now, compute the effective design matrix for the model
    # Multiply X (m x n^2) by X_smooth (n^2 x s), where s is the number of smooth basis functions
    # This gives us an m x s matrix
    
    X_effective <- newX %*% X_smooth
    
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
    data$y = yboot
    # Fit the GAM
    
    fit <- gam(yboot~-1+X_effective, data = data, family = binomial(),
               paraPen = paraPen_list)
    
    # Summary of the model
    summary(fit)
    
    # Extract the estimated image b
    # Compute b = X_smooth %*% coefficients
    b_estimated <- X_smooth %*% coef(fit)
    
    # Reshape b_estimated back to image form
    #estimated_image <- matrix(b_estimated, nx, ny)
    #b_estimated
    #results[[iboot]] = estimated_image
  }
  
  propspline = c()
  propsplinevect = c()
  for(iboot in 1:nboot){
    b_estimated = results[[iboot]]
    propspline = cbind(propspline, b_estimated)
    
    b_estimated = b_estimated[X_matrix!=0]
    propsplinevect = cbind(propsplinevect, b_estimated)
  }
  
  writeMat(spline_filename, propspline=propspline, propsplinevect=propsplinevect)
}
stopCluster(cl)
