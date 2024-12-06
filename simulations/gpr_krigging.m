function smoothed_prop = gpr_krigging(prop, sigma)% Example matrix (propnaive) representing disease distribution

[nrows, ncols] = size(prop);

% Step 1: Logit transformation
logit_transform = @(p) log(p ./ (1 - p));  % Define the logit function
logit_propnaive = logit_transform(prop);  % Apply the logit transformation

% Step 2: Define the spatial grid (coordinates for GP input)
[X1, X2] = meshgrid(1:ncols, 1:nrows);  % Create a grid of coordinates
coords = [X1(:), X2(:)];   % Reshape the grid into a list of coordinate pairs

valuesr = logit_propnaive(:);  % Reshape the logit-transformed matrix into a vector

t_coords = coords(~isnan(prop),:);
t_valuesr = valuesr(~isnan(valuesr));

% Step 3: Fit a Gaussian Process model to the logit-transformed data
gprMdl = fitrgp(t_coords, t_valuesr, 'KernelFunction', 'squaredexponential', ...
     'FitMethod', 'exact', 'PredictMethod', 'exact');

% Step 4: Predict the smoothed logit-transformed values using the Gaussian Process
t_smoothed_logit_values = predict(gprMdl, t_coords);
smoothed_logit_values = NaN(size(coords,1),1);
smoothed_logit_values(~isnan(prop)) = t_smoothed_logit_values;

% Reshape the smoothed values back into a matrix
smoothed_logit_prop = reshape(smoothed_logit_values, nrows, ncols);

% Step 5: Apply the inverse logit transformation to get back to the proportion scale
inv_logit_transform = @(x) exp(x) ./ (1 + exp(x));  % Define the inverse logit function
smoothed_prop = inv_logit_transform(smoothed_logit_prop);  % Apply the inverse logit

