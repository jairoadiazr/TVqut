function smoothed_prop = gpr_krigging_vector(prop, coords)% Example matrix (propnaive) representing disease distribution

logit_transform = @(p) log(p ./ (1 - p));  % Define the logit function
logit_prop = logit_transform(prop);  % Apply the logit transformation

% Step 3: Fit a Gaussian Process model to the logit-transformed data
gprMdl = fitrgp(coords, logit_prop, 'KernelFunction', 'squaredexponential', ...
     'FitMethod', 'exact', 'PredictMethod', 'exact');

% Step 4: Predict the smoothed logit-transformed values using the Gaussian Process
smoothed_logit_prop = predict(gprMdl, coords);

% Step 5: Apply the inverse logit transformation to get back to the proportion scale
inv_logit_transform = @(x) exp(x) ./ (1 + exp(x));  % Define the inverse logit function
smoothed_prop = inv_logit_transform(smoothed_logit_prop);  % Apply the inverse logit

