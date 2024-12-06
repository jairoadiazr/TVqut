
Ns = [500, 5000];
NXs = [30, 50];
TNs = [30, 30*24*4];

snr = 1

niter = 100;

tex_str = '';

bestqut = [];
bestsmooth = [];
bestnaive = [];
bestspline = [];
bestglm = [];


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
        naive_str = '&&empirical';
        glm_str = '&&GLM';
        smooth_str = '&&GPR';
        spline_str = '&&splines';

        rowbestqut = [];
        rowbestnaive = [];
        rowbestspline = [];
        rowbestsmooth = [];
        rowbestglm = [];

        for type = 1:3
            
            if type>1
                qut_str = strcat(qut_str,'&');
                naive_str = strcat(naive_str,'&');
                glm_str = strcat(glm_str,'&');
                smooth_str = strcat(smooth_str,'&');
                spline_str = strcat(spline_str,'&');
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
                propsmoothMSE = [];
                propsplineMSE = [];

                for iter=1:niter

                    filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n_new),'_snr',num2str(snr),'_nx',num2str(nx_new),'_ny',num2str(ny_new),'_tn',num2str(tn_new),'iter',num2str(iter),'.mat')
                    load(filename)
                    

                    %if ~exist('Propspline')
        
                        %Propspline = [];

                        %for num_knots = [4,8]
                        %    for degree = 0:1
                        %        for lambda_spline = [1e-4, 1]
                        %            b_est = binomial_inverse_spline(X, y, nx_new, degree, num_knots, lambda_spline);
                        %            Propspline = [Propspline b_est(:)];
                        %        end
                        %    end
                        %end
                        %save(filename, 'Propspline', '-append');

                    %elseif size(Propspline,2)>8

                        %Propspline = [];

                        %for num_knots = [4,8]
                        %    for degree = 0:1
                        %        for lambda_spline = [1e-4, 1]
                        %            b_est = binomial_inverse_spline(X, y, nx_new, degree, num_knots, lambda_spline);
                        %            Propspline = [Propspline b_est(:)];
                        %        end
                        %    end
                        %end
                        %save(filename, 'Propspline', '-append');

                    %end
                    
                    %for i = 1:size(Propspline,2)
                    %    propspline = Propspline(:,i);
                        Pspline = propspline-min(propspline(:));
                        if max(Pspline(:))~=0
                            Pspline = Pspline/max(Pspline(:));
                        end
                        propsplineMSE = [propsplineMSE, (mean((Pspline(:)-prop(:)).^2,'omitnan'))^0.5];

                    %end
                    %propsplineMSE = [propsplineMSE; propspline_mse];

                    %clear Propspline

                    Phat = prophat-min(prophat(:));
                    if max(Phat(:))~=0
                        Phat = Phat/max(Phat(:));
                    end
                    prophatMSE = [prophatMSE , (mean((Phat(:)-prop(:)).^2))^0.5];
                    
                    Pnaive = propnaive-min(propnaive(:));
                    Pnaive = Pnaive/max(Pnaive(:));
                    propnaiveMSE = [propnaiveMSE , (mean((Pnaive(:)-prop(:)).^2,'omitnan'))^0.5];

                    Psmooth = propsmooth-min(propsmooth(:));
                    Psmooth = Psmooth/max(Psmooth(:));
                    propsmoothMSE = [propsmoothMSE , (mean((Psmooth(:)-prop(:)).^2,'omitnan'))^0.5];

                    if size(X,1)>size(X,2)
                        Pglm = propglm-min(propglm(:));
                        Pglm = Pglm/max(Pglm(:));
                        propglmMSE = [propglmMSE , (mean((Pglm(:)-prop(:)).^2,'omitnan'))^0.5];
                    end
                    
                    % size_prop = size(prop);
                    % A = [prophat(:), ones(length(prophat(:)),1)];
                    % Y = prop;
                    % coeffs = A \ Y(:);
                    % prophat = reshape(A*coeffs, size_prop);
                    % prophatMSE = [prophatMSE , (mean((prophat(:)-prop(:)).^2))^0.5];
                    % 
                    % A = [propsmooth(:), ones(length(propsmooth(:)),1)];
                    % Y = prop;
                    % coeffs = A \ Y(:);
                    % propsmooth = reshape(A*coeffs, size_prop);
                    % propsmoothMSE = [propsmoothMSE , (mean((propsmooth(:)-prop(:)).^2))^0.5];
                    % 
                    % A = [propnaive(~isnan(propnaive)), ones(length(propnaive(~isnan(propnaive))),1)];
                    % Y = prop(~isnan(propnaive));
                    % coeffs = A \ Y(:);
                    % propnaive = A*coeffs;
                    % propnaiveMSE = [propnaiveMSE , (mean((propnaive-prop(~isnan(propnaive))).^2,'omitnan'))^0.5];
                    % 
                    % if size(X,1)>size(X,2)
                    %     A = [propglm(:), ones(length(propglm(:)),1)];
                    %     Y = prop;
                    %     coeffs = A \ Y(:);
                    %     propglm = reshape(A*coeffs, size_prop);
                    %     propglmMSE = [propglmMSE , (mean((propglm(:)-prop(:)).^2))^0.5];
                    % end

                    clear Phat
                    clear prophat
                    clear propglm
                    clear propsmooth
                end
                
                

                prophatSD = round(std(prophatMSE), 2);
                propnaiveSD = round(std(propnaiveMSE), 2);
                propsmoothSD = round(std(propsmoothMSE), 2);
                propsplineSD = round(std(propsplineMSE), 2);

                %[sortprophatMSE, iprophat] = sort(prophatMSE);
                %[sortpronaiveMSE, ipropnaive] = sort(propnaiveMSE);

                %prophatMSE
                
                %prophatMSE = mean(prophatMSE);
                
                [prophatMSE, i_prophat] = sort(prophatMSE);
                [propnaiveMSE, i_propnaive] = sort(propnaiveMSE);
                [propsmoothMSE, i_propsmooth] = sort(propsmoothMSE);
                [propsplineMSE, i_propspline] = sort(propsplineMSE);
                [propglmMSE, i_propglm] = sort(propglmMSE);

                i_bestqut = i_prophat(round(length(i_prophat)/2));
                i_bestnaive = i_propnaive(round(length(i_propnaive)/2));
                i_bestsmooth = i_propsmooth(round(length(i_propsmooth)/2));
                i_bestspline = i_propspline(round(length(i_propspline)/2));

                if length(i_propglm)>0
                    i_bestglm = i_propglm(round(length(i_propglm)/2));
                else
                    i_bestglm = 0;
                end

                prophatMSE = round(mean(prophatMSE), 3);
                propnaiveMSE = round(mean(propnaiveMSE), 3);
                propsmoothMSE = round(mean(propsmoothMSE), 3);
                propsplineMSE = round(mean(propsplineMSE), 3);

                %[propsplineMSE, iminspline] = min(propsplineMSE);
                %iminspline;
                
                if size(X,1)>size(X,2)
                    propglmSD = round(std(propnaiveMSE), 2);
                    propglmMSE = round(mean(propglmMSE), 3);

                    allmse = [prophatMSE, propnaiveMSE, propsmoothMSE, propglmMSE, propsplineMSE];
                    
                    if propglmMSE == min(allmse)
                        glm_str = strcat(glm_str,['&{\bf ' num2str(propglmMSE)], '}');
                    else
                        glm_str = strcat(glm_str,'&',num2str(propglmMSE));
                    end
                else
                    glm_str = strcat(glm_str,'&');
                    allmse = [prophatMSE, propnaiveMSE, propsmoothMSE, Inf, propsplineMSE];
                end
                
                if prophatMSE == min(allmse)
                    qut_str = strcat(qut_str,['&{\bf ' num2str(prophatMSE)], '}');
                else
                    qut_str = strcat(qut_str,'&',num2str(prophatMSE));
                end
                
                if propnaiveMSE == min(allmse)
                    naive_str = strcat(naive_str,['&{\bf ' num2str(propnaiveMSE)], '}');
                else
                    naive_str = strcat(naive_str,'&',num2str(propnaiveMSE));
                end
                
                if propsmoothMSE == min(allmse)
                    smooth_str = strcat(smooth_str,['&{\bf ' num2str(propsmoothMSE)], '}');
                else
                    smooth_str = strcat(smooth_str,'&',num2str(propsmoothMSE));
                end

                if propsplineMSE == min(allmse)
                    spline_str = strcat(spline_str,['&{\bf ' num2str(propsplineMSE)], '}');
                else
                    spline_str = strcat(spline_str,'&',num2str(propsplineMSE));
                end
                
                rowbestqut = [rowbestqut i_bestqut];
                rowbestglm = [rowbestglm i_bestglm];
                rowbestnaive = [rowbestnaive i_bestnaive];
                rowbestspline = [rowbestspline i_bestspline];
                rowbestsmooth = [rowbestsmooth i_bestsmooth];

            end

        end
        qut_str = strcat(qut_str,'\\');
        naive_str = strcat(naive_str,'\\');
        smooth_str = strcat(smooth_str,'\\');
        glm_str = strcat(glm_str,'\\');
        spline_str = strcat(spline_str,'\\');

        bestqut = [bestqut; rowbestqut];
        bestglm = [bestglm; rowbestglm];
        bestnaive = [bestnaive; rowbestnaive];
        bestspline = [bestspline; rowbestspline];
        bestsmooth = [bestsmooth; rowbestsmooth];

        if it==length(TNs)
            tex_str = sprintf('%s\n',tex_str,qut_str, spline_str, smooth_str, naive_str, '\hline');
            
        else
            tex_str = sprintf('%s\n',tex_str,qut_str, spline_str, smooth_str, naive_str);
        end
    end

end    

h=figure;

% subplot(5,6,1)
% filename = strcat('sim_lake_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestqut(1,1)),'.mat');
% load(filename)
% imagesc(prophat);
% set(gca,'XTick',[], 'YTick', [])
% subtitle('(n=500, N=30, t=96)')
% ylabel('TV')
% 
% subplot(5,6,2)
% filename = strcat('sim_lake_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestqut(4,2)),'.mat');
% load(filename)
% imagesc(prophat);
% set(gca,'XTick',[], 'YTick', [])
% subtitle('(n=5000, N=50, t=1)')

subplot(4,4,1)
filename = strcat('sim_river_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestqut(1,3)),'.mat');
load(filename)
imagesc(prophat);
set(gca,'XTick',[], 'YTick', [])
subtitle('(n=500, N=30, t=96)')
ylabel('TV')

subplot(4,4,2)
filename = strcat('sim_river_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestqut(4,4)),'.mat');
load(filename)
imagesc(prophat);
set(gca,'XTick',[], 'YTick', [])
subtitle('(n=5000, N=50, t=1)')

subplot(4,4,3)
filename = strcat('sim_sides_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestqut(1,5)),'.mat');
load(filename)
imagesc(prophat);
set(gca,'XTick',[], 'YTick', [])
subtitle('(n=500, N=30, t=96)')

subplot(4,4,4)
filename = strcat('sim_sides_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestqut(4,6)),'.mat');
load(filename)
imagesc(prophat);
set(gca,'XTick',[], 'YTick', [])
subtitle('(n=5000, N=50, t=1)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(5,6,7)
% filename = strcat('sim_lake_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestnaive(1,1)),'.mat');
% load(filename)
% imagesc(propnaive);
% set(gca,'XTick',[], 'YTick', [])
% ylabel('empirical')
% 
% subplot(5,6,8)
% filename = strcat('sim_lake_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestnaive(4,2)),'.mat');
% load(filename)
% imagesc(propnaive);
% set(gca,'XTick',[], 'YTick', [])

subplot(4,4,5)
filename = strcat('sim_river_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestnaive(1,3)),'.mat');
load(filename)
imagesc(propspline);
set(gca,'XTick',[], 'YTick', [])
ylabel('splines')

subplot(4,4,6)
filename = strcat('sim_river_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestnaive(4,4)),'.mat');
load(filename)
imagesc(propspline);
set(gca,'XTick',[], 'YTick', [])

subplot(4,4,7)
filename = strcat('sim_sides_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestnaive(1,5)),'.mat');
load(filename)
imagesc(propspline);
set(gca,'XTick',[], 'YTick', [])

subplot(4,4,8)
filename = strcat('sim_sides_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestnaive(4,6)),'.mat');
load(filename)
imagesc(propspline);
set(gca,'XTick',[], 'YTick', [])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(5,6,13)
% filename = strcat('sim_lake_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestsmooth(1,1)),'.mat');
% load(filename)
% imagesc(propsmooth);
% set(gca,'XTick',[], 'YTick', [])
% ylabel('GPR')
% 
% subplot(5,6,14)
% filename = strcat('sim_lake_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestsmooth(4,2)),'.mat');
% load(filename)
% imagesc(propsmooth);
% set(gca,'XTick',[], 'YTick', [])

subplot(4,4,9)
filename = strcat('sim_river_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestsmooth(1,3)),'.mat');
load(filename)
imagesc(propsmooth);
set(gca,'XTick',[], 'YTick', [])
ylabel('GPR')

subplot(4,4,10)
filename = strcat('sim_river_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestsmooth(4,4)),'.mat');
load(filename)
imagesc(propsmooth);
set(gca,'XTick',[], 'YTick', [])

subplot(4,4,11)
filename = strcat('sim_sides_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestsmooth(1,5)),'.mat');
load(filename)
imagesc(propsmooth);
set(gca,'XTick',[], 'YTick', [])

subplot(4,4,12)
filename = strcat('sim_sides_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestsmooth(4,6)),'.mat');
load(filename)
imagesc(propsmooth);
set(gca,'XTick',[], 'YTick', [])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(5,6,19)
% filename = strcat('sim_lake_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestspline(1,1)),'.mat');
% load(filename)
% imagesc(propspline);
% set(gca,'XTick',[], 'YTick', [])
% ylabel('splines')
% 
% subplot(5,6,20)
% filename = strcat('sim_lake_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestspline(4,2)),'.mat');
% load(filename)
% imagesc(propspline);
% set(gca,'XTick',[], 'YTick', [])

subplot(4,4,13)
filename = strcat('sim_river_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestspline(1,3)),'.mat');
load(filename)
imagesc(propnaive);
set(gca,'XTick',[], 'YTick', [])
ylabel('empirical')

subplot(4,4,14)
filename = strcat('sim_river_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestspline(4,4)),'.mat');
load(filename)
imagesc(propnaive);
set(gca,'XTick',[], 'YTick', [])

subplot(4,4,15)
filename = strcat('sim_sides_bird_moves_n500_snr1_nx30_ny30_tn30','iter',num2str(bestspline(1,5)),'.mat');
load(filename)
imagesc(propnaive);
set(gca,'XTick',[], 'YTick', [])

subplot(4,4,16)
filename = strcat('sim_sides_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestspline(4,6)),'.mat');
load(filename)
imagesc(propnaive);
set(gca,'XTick',[], 'YTick', [])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% empty_image = zeros(30,30);
% % subplot(5,6,25)
% % imagesc(empty_image);
% % set(gca,'XTick',[], 'YTick', [])
% % ylabel('GLM')
% % 
% % subplot(5,6,26)
% % filename = strcat('sim_lake_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestglm(4,2)),'.mat');
% % load(filename)
% % imagesc(propglm);
% % set(gca,'XTick',[], 'YTick', [])
% 
% subplot(5,4,17)
% imagesc(empty_image);
% set(gca,'XTick',[], 'YTick', [])
% ylabel('GLM')
% 
% subplot(5,4,18)
% filename = strcat('sim_river_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestglm(4,4)),'.mat');
% load(filename)
% imagesc(propglm);
% set(gca,'XTick',[], 'YTick', [])
% 
% subplot(5,4,19)
% imagesc(empty_image);
% set(gca,'XTick',[], 'YTick', [])
% 
% subplot(5,4,20)
% filename = strcat('sim_sides_bird_moves_n5000_snr1_nx50_ny50_tn2880','iter',num2str(bestglm(4,6)),'.mat');
% load(filename)
% imagesc(propglm);
% set(gca,'XTick',[], 'YTick', [])


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'position',[pos(1),pos(2),pos(3)*1.3,pos(4)*1.7])
annotation('line', [0.52 0.52], [0.1 .95], 'Color', 'k', 'LineWidth', 1.5);
annotation('textbox', [0.22, 0.93, 0.18, 0.05], 'String', 'river', 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 11);
annotation('textbox', [0.635, 0.93, 0.18, 0.05], 'String', 'lake+corner', 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 11);

pos = get(h,'Position');
set(h, 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(h,'PaperPositionMode','Auto')
print(h,'asymptotics','-dpdf','-r0','-bestfit')