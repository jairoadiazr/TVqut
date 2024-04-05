function [lambda_qut, lambda_max] = lambdaQUT(X, y, D_TV, MC)
    y01=y;
    X0 = X*ones(size(X,2),1);
    x0 = glmfit(X0,y,'binomial','Constant','off');

    mu0 = x0*X0;
    p0=1./(1+exp(-mu0));

    parfor i = 1:MC
        i
        y0_temp = binornd(1,p0);
        [w, lambda] = SolveInfinityNormDeer_cvx(y0_temp,X,D_TV);
        Lambdas(i) = lambda;
    end
    lambda_qut = quantile(Lambdas, 0.95);
    lambda = lambda_qut;
    
    [w,lambda_max] = SolveInfinityNormDeer_cvx(y01,X,D_TV);