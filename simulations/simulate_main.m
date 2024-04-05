
%simulate the main movements
%type = 1; %lake
%type = 2; %river
%type = 3; %sides

flag = 0;
niter = 100;

for type = 1:3

    if type == 1
        strtype = 'lake'
    elseif type ==2
        strtype = 'river'
    elseif type ==3
        strtype = 'sides'
    end
    
    for iter = 1:niter
        %SNR
        snr = 1;
        
        filename = strcat(strtype,'_bird_moves_main','_iter',num2str(iter),'.mat');
        if ~isfile(filename)
            simulate_movements_main
        end
         
        %simulate X and y
        filename = strcat(strtype,'_bird_moves_main_snr',num2str(snr),'iter',num2str(iter),'.mat');
        if ~isfile(filename)
            simulate_X_y_main
        end

        %parameter sets
        NXs = [10, 30, 50];
        Ns = [500, 2000, 5000];
        TNs = [6, 30, 30*24*4];
                
        for nx_new = NXs
        
            for n_new = Ns
                
                for tn_new = TNs
                    ny_new = nx_new;
                    
                    filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n_new),'_snr',num2str(snr),'_nx',num2str(nx_new),'_ny',num2str(ny_new),'_tn',num2str(tn_new),'iter',num2str(iter),'.mat')
                    if isfile(filename)
                        
                        continue
                    end

                    [n_new, nx_new, tn_new]
                    
                    'reading main file'
                    filename = strcat(strtype,'_bird_moves_main_snr',num2str(snr),'iter',num2str(iter),'.mat');
                    load(filename)
        
                    %reshape movements from main
                    simulate_reshape_movements_from_main

                    %calculate tv and naive estimator
                    simulate_estimate_prop
                    
                    filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n_new),'_snr',num2str(snr),'_nx',num2str(nx_new),'_ny',num2str(ny_new),'_tn',num2str(tn_new),'iter',num2str(iter),'.mat')
                    load(filename)

                    %calculate GLM estimator
                    if size(X,1) > size(X,2)
                        propglm = glmfit(X/tn_new,y,'binomial','Constant','off');
                        propglm = reshape(propglm,nx_new,ny_new);
                    else
                        propglm = [];
                    end
                    save(filename,'propglm','-append')
                end
            end
        end
    end
end