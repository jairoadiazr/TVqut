
maxiter = 100;
weight = 1000;

interpolations = [15];
OVs = [0];
year1s = [0];
ehdvs = [1, 2, 6];

calculate_all = 1;


for iter = 1:maxiter
    '#####################'
    iter

    rng(iter)
    if calculate_all == 1
    
        for interpolation = interpolations
        
            for OV = OVs
        
                if OV == 0
                    ov = '';
                else
                    ov = 'OV';
                end
        
                for year1 = year1s
                
                    if year1 == 0
                        year = 'ALL';
                    else
                        year = year1;
                    end
            
                    for ehdv = ehdvs
                                
                        y=csvread(strcat('Y_vector_EHDV-',num2str(ehdv),'_weight_',num2str(weight),'_year_',num2str(year), '_interpolation_',num2str(interpolation), '_', ov, '.csv'), 1, 1);
                        X=csvread(strcat('X_matrix_EHDV-',num2str(ehdv),'_weight_',num2str(weight),'_year_',num2str(year), '_interpolation_',num2str(interpolation), '_', ov, '.csv'), 1, 1);
            
                        %if cap == 1
                        %    X = X.*(X<repmat(floor(quantile(X',0.90))',1,size(X,2))) + repmat(floor(quantile(X',0.90))',1,size(X,2)).*(X>=repmat(floor(quantile(X',0.90))',1,size(X,2)));
                        %end
                
                        D=csvread(strcat('D_matrix_EHDV-',num2str(ehdv),'_weight_',num2str(weight),'_year_',num2str(year), '_interpolation_',num2str(interpolation), '_', ov, '.csv'), 1, 1);
                        
                        y=y(sum(X')>1000);
                        X=X(sum(X')>1000,:);
        
                        [X, y] = bootstrap_moves_plus(X, y, 600, 2000);
                        
                        D=csvread(strcat('D_matrix_EHDV-',num2str(ehdv),'_weight_',num2str(weight),'_year_',num2str(year), '_interpolation_',num2str(interpolation), '_', ov, '.csv'), 1, 1);
                        points = csvread(strcat('points_EHDV-',num2str(ehdv),'_weight_',num2str(weight),'_year_',num2str(year), '_interpolation_',num2str(interpolation), '_', ov, '.csv'), 1, 1);
                        points = round(points*weight-min(points*weight)+1);
                        
			    
					    %CALCULATE QUT estimator
					    
                        D_TV=D;
                        D_TV=sparse(D_TV);
                        
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
                        

                        prop_vect_hat = sol.x;
                        prop_hat = NaN(max(points));
                        for i = 1:length(sol.x)
                            prop_hat(points(i,1), points(i,2)) = prop_vect_hat(i);
                        end
                        
                        prophat = sol.x;
                        
                        propnaive = sum(X(y==1,:))./sum(X);
                        prop_vect_naive = propnaive;
                        prop_naive = NaN(max(points));
                        for i = 1:length(sol.x)
                            prop_naive(points(i,1), points(i,2)) = prop_vect_naive(i);
                        end
                        
                        propglm = glmfit(X,y,'binomial','constant','off');
                        prop_vect_glm = propglm;
                        prop_glm = NaN(max(points));
                        for i = 1:length(sol.x)
                            prop_glm(points(i,1), points(i,2)) = prop_vect_glm(i);
                        end
    
                        propsmooth = gpr_krigging_vector(propnaive, points);
                        prop_vect_smooth = propsmooth;
                        prop_smooth = NaN(max(points));
                        for i = 1:length(propnaive)
                            prop_smooth(points(i,1), points(i,2)) = prop_vect_smooth(i);
                        end
    
                        filename = strcat('deer_EHDV-',num2str(ehdv),'_weight_',num2str(weight),'_year_',num2str(year), '_interpolation_',num2str(interpolation), '_', ov,'_PLUS_iter', num2str(iter),'.mat')
                        save(filename,'y','prop_vect_naive','prop_vect_glm','prop_vect_hat','prop_vect_smooth','propsmooth','prop_smooth','prophat','prop_hat','propnaive','prop_naive','propglm','prop_glm','lambda_max','lambda_qut','X','D_TV','-v7.3')
                        
                        %points = csvread(strcat('points_EHDV-',num2str(ehdv),'_weight_',num2str(weight),'_year_',num2str(year), '_interpolation_',num2str(interpolation), '_', ov,'.csv'), 1, 1);
                        %points(:,3) = prop_vect_hat;
                        %writematrix(points, strcat('result_coords_EHDV-',num2str(ehdv),'_weight_',num2str(weight),'_year_',num2str(year), '_interpolation_',num2str(interpolation), '_', ov,'.txt'))
                    end
                end
            end
        end
    end
end