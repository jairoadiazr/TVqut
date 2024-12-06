
%simulate the main movements
%type = 1; %lake
%type = 2; %river
%type = 3; %sidesj nvbhjikbmnyuh

flag = 0;

%for type = 1:3
for type = [1]

    
    if type == 1
        strtype = 'lake'
    elseif type ==2
        strtype = 'river'
    elseif type ==3
        strtype = 'sides'
    end
    
    filename = strcat('test_',strtype,'_bird_moves_main.mat')
    if isfile(filename)
        load(filename)
    else
        simulate_movements_main_test
    end


    prop_original = prop;
    nx_original = nx;
    ny_original = ny;
    n_original = n;
    tn_original = tn;
    X_original = X;

    for iter = 1:1000
        
        for snr = [0,1]
        
            propvec_original = reshape(snr*prop_original/100, nx_original*ny_original, 1);
            mu=X_original * propvec_original;
            p_original=1./(1+exp(-mu));
            y_original=binornd(1,p_original);
            propvec_original = reshape(propvec_original,nx_original,ny_original);
            
            %NXs = [10,30,50];
            NXs = [10];
            Ns = 50:50:400;
            
            for nx_new = NXs
                ny_new = nx_new;

                for n_new = Ns
                    tn_new = 30;
                    
                    filename = strcat('test_sim_',strtype,'_bird_moves_n',num2str(n_new),'_snr',num2str(snr),'_nx',num2str(nx_new),'_ny',num2str(ny_new),'_tn',num2str(tn_new),'iter',num2str(iter),'.mat')
        
                    %reshape movements from main
                    tn_new = 30;

                    tn = tn_original;
                    nx = nx_original;
                    ny = ny_original;
                    n = n_original;
                    y = y_original;
                    p = p_original;
                    propvec = propvec_original;
                    
                    simulate_reshape_movements_from_main

                    %%%Perform test using lambdaQUT
                    simulate_test_prop
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %Perform other tests
                    if n_new > (nx*ny)

                        [b,dev] = glmfit(X,y,'binomial','Constant','off');
                        X0 = X*ones(size(X,2),1);
                        [b0,dev0] = glmfit(X0,y,'binomial','Constant','off');
    
                        D = dev0 - dev;
    
                        df = size(X,2) - size(X0,2);
    
                        pvalue_llratio = 1 - chi2cdf(D,df);
                        llratio = pvalue_llratio<0.05;
                        
                        mu0 = b0*X0;
                        p0=1./(1+exp(-mu0));
                        
                        MC = 100;
                        Ds = zeros(MC,1);
    
                        parfor i=1:MC
                            y0_temp = binornd(1,p0);
                            [~,dev_temp] = glmfit(X,y0_temp,'binomial','Constant','off');
                            [~,dev0_temp] = glmfit(X0,y0_temp,'binomial','Constant','off');
                            D_temp = dev0_temp - dev_temp;
                            Ds(i) = D_temp;
                        end
                        D_threshold = quantile(Ds,0.95);
                        llratio_exact = D>D_threshold;
                        save(filename,'pvalue_llratio','llratio', 'llratio_exact', 'Ds', 'D','-append')
                        [llratio, llratio_exact]
                    end
                end
            end
        end
    end
end