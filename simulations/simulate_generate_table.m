
Ns = [500, 2000, 5000];
NXs = [10, 30, 50];
TNs = [6, 30, 30*24*4];

snr = 1

niter = 50;

tex_str = '';
for inx = 1:length(NXs)

    nx_new = NXs(inx);
    ny_new = nx_new;
    
    for it = 1:length(TNs)
        tn_new = TNs(it);
            
        if it == 1
            qut_str = strcat('$',num2str(nx_new),'$&');
        else
            qut_str = '&';
        end
        Tstr = 30*24*4/tn_new;
        qut_str = strcat(qut_str,num2str(Tstr),'&TV');
        naive_str = '&&naive';
        glm_str = '&&GLM';

        for type = 1:3
            
            if type>1
                qut_str = strcat(qut_str,'&');
                naive_str = strcat(naive_str,'&');
                glm_str = strcat(glm_str,'&');
            end
            
            prop=zeros(nx_new,ny_new);

            if type == 1
                strtype = 'lake';
                strtype2 = 'lake';
                for i = 1:nx_new
                    for j = 1:ny_new
                        if (i-4*ny_new/7)^2 + (j-4*ny_new/7)^2 <= (ny_new/5)^2
                            prop(i,j)=1;
                        end
                    end
                end
            
            elseif type == 2
            
                strtype = 'river';
                strtype2 = 'river';
            %prop=zeros(nx,ny);
            
               for i = 1:nx_new
                    for j = 1:ny_new
                        if i>(2*j-3*ny_new/5)&&(i<2*j-ny_new/6)
                            prop(i,j)=1;
                        end
                    end
               end
               %mean(prop(:))
               %imagesc(prop)
            
            elseif type ==3
            
                strtype = 'sides';
                strtype2 = 'lake + corner';
                for i = 1:nx_new
                    for j = 1:ny_new
                        if (i-ny_new/3)^2 + (j-2*ny_new/3)^2 <= (ny_new/5)^2
                            prop(i,j)=1;
                        end
            
                        %if j>(i-2*ny/3)&&(j<i-4*ny/7)
                        if (j<i-5*ny_new/7)
                            prop(i,j)=1;
                        end
                    end
                end
            
            end
            
            for in = 1:length(Ns)

                n_new = Ns(in);

                prophatMSE = [];
                propnaiveMSE = [];
                propglmMSE = [];

                for iter=1:niter
                    filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n_new),'_snr',num2str(snr),'_nx',num2str(nx_new),'_ny',num2str(ny_new),'_tn',num2str(tn_new),'iter',num2str(iter),'.mat');
                    load(filename)

                    Phat = prophat-min(prophat(:));
                    if max(Phat(:))~=0
                        Phat = Phat/max(Phat(:));
                    end
                    prophatMSE = [prophatMSE , (mean((Phat(:)-prop(:)).^2))^0.5];
                    
                    Pnaive = propnaive-min(propnaive(:));
                    Pnaive = Pnaive/max(Pnaive(:));
                    propnaiveMSE = [propnaiveMSE , (mean((Pnaive(:)-prop(:)).^2,'omitnan'))^0.5];

                    
                    if size(X,1)>size(X,2)
                        Pglm = propglm-min(propglm(:));
                        Pglm = Pglm/max(Pglm(:));
                        propglmMSE = [propglmMSE , (mean((Pglm(:)-prop(:)).^2,'omitnan'))^0.5];
                    end
                    
                    clear Phat
                    clear prophat
                    clear propglm
                end
                
                

                prophatSD = round(std(prophatMSE), 2);
                propnaiveSD = round(std(propnaiveMSE), 2);

                %[sortprophatMSE, iprophat] = sort(prophatMSE);
                %[sortpronaiveMSE, ipropnaive] = sort(propnaiveMSE);

                %prophatMSE
                
                %prophatMSE = mean(prophatMSE);
                prophatMSE = round(mean(prophatMSE), 3);
                propnaiveMSE = round(mean(propnaiveMSE), 3);

                if size(X,1)>size(X,2)
                    propglmSD = round(std(propnaiveMSE), 2);
                    propglmMSE = round(mean(propglmMSE), 3);
                    glm_str = strcat(glm_str,'&',num2str(propglmMSE));
                else
                    glm_str = strcat(glm_str,'&');
                end

                qut_str = strcat(qut_str,'&',num2str(prophatMSE));
                naive_str = strcat(naive_str,'&',num2str(propnaiveMSE));

            end

        end
        qut_str = strcat(qut_str,'\\');
        naive_str = strcat(naive_str,'\\');
        glm_str = strcat(glm_str,'\\');

        
        if it==length(TNs)
            tex_str = sprintf('%s\n',tex_str,qut_str,glm_str, naive_str,'\hline')
        else
            tex_str = sprintf('%s\n',tex_str,qut_str, glm_str, naive_str);
        end
    end

end    