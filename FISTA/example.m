load('deer_EHDV-6_weight_1000_year_ALL_interpolation_15_.mat')

X0 = X*ones(size(X,2),1);
options.x0 = x0*ones(size(X,2),1);                
options.maxItr = 10000;
options.fistaitr = 10000;
options.tol_l1 = 1e-4;
options.tol_Q2 = 1e-6;
options.tol_gs = 1e-6;
options.tol_TV = 1e-20;
options.fistatol = 1e-7;

sol = MFISTAlosses_TV(X, y, 'entropy', lambda_qut, D_TV, 2, [],options);

prop_vect_hat = sol.x;
prop_hat = NaN(max(points));
for i = 1:length(sol.x)
    prop_hat(points(i,1), points(i,2)) = prop_vect_hat(i);
end
imagesc(prop_hat)