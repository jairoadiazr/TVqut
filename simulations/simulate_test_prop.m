

filename = strcat('test_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn_new),'iter',num2str(iter),'.mat')

X0 = X*ones(size(X,2),1);
x0 = glmfit(X0,y,'binomial','Constant','off');
propvec0_hat = x0*ones(size(X,2),1);
mu0 = x0*X0;
p0=1./(1+exp(-mu0));

options.x0 = propvec0_hat/tn_new;

%lambdaqut
MC = 100;
Lambdas = zeros(MC,1);
Lambdas_lasso = zeros(MC,1);

%'calculating lambda qut'


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

        X0 = X*ones(size(X,2),1);
        x0_temp = glmfit(X0,y0_temp,'binomial','Constant','off');
        mu0_temp = x0_temp*X0;
        p0_temp = 1./(1+exp(-mu0_temp));
        lambda_lasso = max(abs(X'*(y0_temp - p0_temp)));
        Lambdas_lasso(i) = lambda_lasso;

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
lambda_lasso_qut = quantile(Lambdas_lasso, 0.95);
lambda_lasso_max = max(abs(X'*(y - p0)));

lambda = lambda_qut;

[w,lambda_max] = SolveInfinityNormDeer_cvx(y,X/tn_new,D_TV);

if isnan(lambda_max)
    %lambda_max
    [w,lambda_max] = SolveInfinityNormDeer_cvx(y,X/tn_new/10,D_TV);
    lambda_max = lambda_max * 10;
    %lambda_max
end

if isnan(lambda_max)
    %lambda_max
    [w,lambda_max] = SolveInfinityNormDeer_cvx(y,X/tn_new*10,D_TV);
    lambda_max = lambda_max / 10;
    %lambda_max
end

if isnan(lambda_max)
    afsdfsd
end

[iter, lambda_qut, lambda_max]

propnaive = reshape(sum(X(y==1,:))./sum(X),nx,ny);
filename = strcat('test_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn_new),'iter',num2str(iter),'.mat');
save(filename,'Lambdas','lambda_max','lambda_qut','lambda_lasso_max', 'lambda_lasso_qut', 'Lambdas_lasso','X','D_TV','n','nx','ny','y','tn_new','-v7.3')
