%FISTA OPTIONS
options.maxItr = 20000;
options.fistaitr = 20000;
options.tol_l1 = 1e-4;
options.tol_Q2 = 1e-6;
options.tol_gs = 1e-6;
options.tol_TV = 1e-20;
options.fistatol = 1e-6;
tol_TV=options.tol_TV;
maxItr=options.maxItr;



X0 = X*ones(size(X,2),1);
x0 = glmfit(X0,y,'binomial','Constant','off');
propvec0_hat = x0*ones(size(X,2),1);
mu0 = x0*X0;
p0=1./(1+exp(-mu0));

options.x0 = propvec0_hat/tn_new;

%lambdaqut
MC = 100;
Lambdas = zeros(MC,1);

'calculating lambda qut'


if n<10000
    parfor i = 1:MC
        y0_temp = binornd(1,p0);
        [~, lambda] = SolveInfinityNormDeer_cvx(y0_temp,X/tn_new,D_TV);

        if isnan(lambda)
            i
            lambda
            [~, lambda] = SolveInfinityNormDeer_cvx(y0_temp,X/tn_new/10,D_TV);
            lambda = lambda * 10;
            lambda
        end
        
        if isnan(lambda)
            i
            lambda
            [~, lambda] = SolveInfinityNormDeer_cvx(y0_temp,X/tn_new*10,D_TV);
            lambda = lambda / 10;
            lambda
        end
        Lambdas(i) = lambda;
    end
else
    for i = 1:MC
        i
        y0_temp = binornd(1,p0);
        [~, lambda] = SolveInfinityNormDeer_cvx(y0_temp,X/tn_new,D_TV);

        if isnan(lambda)
            [~, lambda] = SolveInfinityNormDeer_cvx(y0_temp,X/tn_new/10,D_TV);
            lambda = lambda * 10;
        end
        
        if isnan(lambda)
            [~, lambda] = SolveInfinityNormDeer_cvx(y0_temp,X/tn_new*10,D_TV);
            lambda = lambda / 10;
        end

        Lambdas(i) = lambda;
    end
end    

lambda_qut = quantile(Lambdas, 0.95);
lambda = lambda_qut;

[w,lambda_max] = SolveInfinityNormDeer_cvx(y,X/tn_new,D_TV);

if isnan(lambda_max)
    lambda_max
    [w,lambda_max] = SolveInfinityNormDeer_cvx(y,X/tn_new/10,D_TV);
    lambda_max = lambda_max * 10;
    lambda_max
end

if isnan(lambda_max)
    lambda_max
    [w,lambda_max] = SolveInfinityNormDeer_cvx(y,X/tn_new*10,D_TV);
    lambda_max = lambda_max / 10;
    lambda_max
end

if isnan(lambda_max)
    afsdfsd
end

if lambda_qut<lambda_max
    sol = MFISTAlosses_TV(X/tn_new, y, 'entropy', lambda, D_TV, 2, [], options);
    prophat=reshape(sol.x,nx,ny);
else
    prophat = x0*ones(nx,ny);
end
[lambda_qut, lambda_max]



propnaive = reshape(sum(X(y==1,:))./sum(X),nx,ny);
filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn_new),'iter',num2str(iter),'.mat')
save(filename,'prophat','propnaive','Lambdas','lambda_max','lambda_qut','X','D_TV','n','nx','ny','y','propvec','tn','-v7.3')
