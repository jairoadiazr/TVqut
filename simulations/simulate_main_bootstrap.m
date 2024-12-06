NXs = [30];
Ns = [500];
maxiter = 20;
nboot = 50;
tn = 2880;
snr = 1;


for iter = 1:maxiter
    
    for type = [1:3]

        if type == 1
            strtype = 'lake'
        elseif type ==2
            strtype = 'river'
        elseif type ==3
            strtype = 'sides'
        end

        for n = Ns
            for nx = NXs
                ny = nx;
                
                
                for iboot = 1:nboot 
                    iboot

                    filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(iter),'.mat');
                    load(filename)

                    rng(iboot)
                    n_replicas = 2*n;
                    n_movements = 1000;

                    [X, y] = bootstrap_moves_plus(X, y, n_replicas, n_movements);

                    X0 = X*ones(size(X,2),1);
                    x0 = glmfit(X0,y,'binomial','Constant','off');
                    options.x0 = x0*ones(size(X,2),1);
                    
                    propvec0_hat = x0*ones(size(X,2),1);
                    mu0 = x0*X0;
                    p0=1./(1+exp(-mu0));
                    

				    %FISTA OPTIONS
                    options.maxItr = 100000;
                    options.fistaitr = 100000;
                    
                    options.tol_l1 = 1e-5;
                    options.tol_Q2 = 1e-7;
                    options.tol_gs = 1e-7;
                    options.tol_TV = 1e-21;
                    
                    tol_TV=options.tol_TV;
                    maxItr=options.maxItr;
                    %x0_tv=zeros(nx*ny,1);
                    
                    options.fistatol = 1e-8;
                    
                    MC = 100;

                    %Calculate TV QUT estimator
                    [sol, lambda_qut, lambda_max] = tvQUT(X,y,D_TV,MC,options);                    
                    

                    prophat_boot=reshape(sol.x,nx,ny);
                    propnaive_boot = reshape(sum(X(y==1,:))./sum(X),nx,ny);
                    propsmooth_boot = gpr_krigging(propnaive_boot);

                    propglm_boot = glmfit(X,y,'binomial','constant','off');
                    propglm_boot = reshape(propglm_boot,nx,ny);
                    lambda_qut_boot = lambda_qut;

                    filename = strcat('BOOT_', num2str(iboot), + 'sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(iter),'.mat');
                    save(filename,'lambda_qut_boot', 'prophat_boot','propnaive_boot','propsmooth_boot','propglm_boot')
                end
            end
        end
    end
end