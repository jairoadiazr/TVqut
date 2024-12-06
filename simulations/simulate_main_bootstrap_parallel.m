NXs = [30];
Ns = [500];
maxiter = 50;
nboot = 48;
tn = 2880;
snr = 1;

for iter = 1:maxiter
    

    for type = [1:3]
        strcat('iter',num2str(iter))
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
                n_replicas = 5000;
                n_movements = 720;

                save_filename = strcat('BOOT_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(iter),'nboot',num2str(nboot),'tn_new', num2str(n_movements), 'n_new', num2str(n_replicas),'.mat');
                
                if exist(save_filename, 'file') == 2
                    continue
                else
                    save_filename
                end

                filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(iter),'.mat');
                load(filename)



			    %FISTA OPTIONS
                simoptions.maxItr = 100000;
                simoptions.fistaitr = 100000;
                
                simoptions.tol_l1 = 1e-5;
                simoptions.tol_Q2 = 1e-7;
                simoptions.tol_gs = 1e-7;
                simoptions.tol_TV = 1e-21;
                
                tol_TV=simoptions.tol_TV;
                maxItr=simoptions.maxItr;
                %x0_tv=zeros(nx*ny,1);
                
                simoptions.fistatol = 1e-8;
                
                MC = 100
                
                all_prophat_boot = zeros(nx*ny, nboot);
                all_propsmooth_boot = zeros(nx*ny, nboot);
                all_propnaive_boot = zeros(nx*ny, nboot);
                %all_propglm_boot = zeros(nx*ny, nboot);
                all_lambdaqut_boot = zeros(1,nboot);
                parfor iboot = 1:nboot 
                    disp(['Completed boot: ', num2str(iboot)]);

                    rng(iboot)
                    options = simoptions;

                    [Xboot, yboot] = bootstrap_moves_plus(X, y, n_replicas, n_movements);

                    
                    options.x0 = Xboot*ones(size(Xboot,2),1);
                    options.x0 = glmfit(options.x0,yboot,'binomial','Constant','off');
                    options.x0 = options.x0*ones(size(Xboot,2),1);
                    
                    %Calculate TV QUT estimator
                    [prophat_boot, lambda_qut, lambda_max] = tvQUT(Xboot,yboot,D_TV,MC,options);                    
                    

                    prophat_boot=reshape(prophat_boot.x,nx,ny);
                    propnaive_boot = reshape(sum(Xboot(yboot==1,:))./sum(Xboot),nx,ny);
                    propsmooth_boot = gpr_krigging(propnaive_boot);

                    %propglm_boot = glmfit(Xboot,yboot,'binomial','constant','off');
                    %propglm_boot = reshape(propglm_boot,nx,ny);
                    lambda_qut_boot = lambda_qut;
                    
                    all_lambdaqut_boot(iboot) = lambda_qut_boot

                    all_prophat_boot(:,iboot) = prophat_boot(:);
                    all_propsmooth_boot(:,iboot) = propsmooth_boot(:);
                    all_propnaive_boot(:,iboot) = propnaive_boot(:);
                    %all_propglm_boot(:,iboot) = propglm_boot(:);
                
                    
                end
                filename = strcat('BOOT_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(iter),'nboot',num2str(nboot),'tn_new', num2str(n_movements), 'n_new', num2str(n_replicas),'.mat');
                %save(filename,'all_prophat_boot', 'all_propsmooth_boot','all_propnaive_boot','all_propglm_boot','all_lambdaqut_boot')
                save(filename,'all_prophat_boot', 'all_propsmooth_boot','all_propnaive_boot','all_lambdaqut_boot')
                clear all_prophat_boot
                clear all_propsmooth_boot
                clear all_propnaive_boot
                %clear all_propglm_boot
            end
        end
    end
end


iboot = 2
rng(iboot)
options = simoptions;

[Xboot, yboot] = bootstrap_moves_plus(X, y, n_replicas, n_movements);
propglm_boot = glmfit(Xboot,yboot,'binomial','constant','off');