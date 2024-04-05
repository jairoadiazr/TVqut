NXs = [10, 30, 50];
h = figure

n_new = 2000;
snr = 1;


niter = 100;

strtype = 'sides';
strtype2 = 'lake + corner'

for inx = 1:length(NXs)
    nx_new = NXs(inx);
    tn_new = 30;

    prop=zeros(nx_new,nx_new);
    for i = 1:nx_new
        for j = 1:nx_new
            if (i-nx_new/3)^2 + (j-2*nx_new/3)^2 <= (nx_new/5)^2
                prop(i,j)=1;
            end
    
            %if j>(i-2*ny/3)&&(j<i-4*ny/7)
            if (j<i-5*nx_new/7)
                prop(i,j)=1;
            end
        end
    end

    prophatMSE = [];
    propnaiveMSE = [];
    propglmMSE = [];

    for iter=1:niter
        filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n_new),'_snr',num2str(snr),'_nx',num2str(nx_new),'_ny',num2str(nx_new),'_tn',num2str(tn_new),'iter',num2str(iter),'.mat');
        load(filename)
    
        Phat = prophat-min(prophat(:));
        if max(Phat(:))~=0
            Phat = Phat/max(Phat(:));
        end
        Phat = Phat/max(Phat(:));
        prophatMSE = [prophatMSE , (mean((Phat(:)-prop(:)).^2))^0.5];
        
        if size(X,1) > size(X,2)
            Pglm = propglm-min(propglm(:));
            if max(Pglm(:))~=0
                Pglm = Pglm/max(Pglm(:));
            end
            Pglm = Pglm/max(Pglm(:));
            propglmMSE = [propglmMSE , (mean((Pglm(:)-prop(:)).^2))^0.5];
        end
    
        Pnaive = propnaive-min(propnaive(:));
        Pnaive = Pnaive/max(Pnaive(:));
        propnaiveMSE = [propnaiveMSE , (mean((Pnaive(:)-prop(:)).^2,'omitnan'))^0.5];
    
    end

    [sortprophatMSE, iprophat] = sort(prophatMSE);
    [sortpronaiveMSE, ipropnaive] = sort(propnaiveMSE);
    
    if size(X,1)>size(X,2)
        [sortpropglmMSE, ipropglm] = sort(propglmMSE);
    end

    filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n_new),'_snr',num2str(snr),'_nx',num2str(nx_new),'_ny',num2str(nx_new),'_tn',num2str(tn_new),'iter',num2str(iprophat(niter/2)),'.mat')
    load(filename)
    subplot(3, 3, inx)          
    imagesc(prophat)
    set(gca,'XTick',[], 'YTick', [])
    %xlabel(strcat('median MSE=',num2str(round(prophatMSE(iprophat(niter/2)),3))))
    title(strcat('N=',num2str(nx)))
    
    if inx==1
        ylabel('TV')
    end

    if size(X,1) > size(X,2)
        filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n_new),'_snr',num2str(snr),'_nx',num2str(nx_new),'_ny',num2str(nx_new),'_tn',num2str(tn_new),'iter',num2str(ipropglm(niter/2)),'.mat')
        load(filename)
        subplot(3, 3, 3 + inx)          
        imagesc(propglm)
        set(gca,'XTick',[], 'YTick', [])
        %xlabel(strcat('median MSE=',num2str(round(propglmMSE(ipropglm(niter/2)),3))))
    else
        subplot(3, 3, 3 + inx) 
        imshow(prop.*0)
    end

    if inx==1
        ylabel('GLM')
    end

    filename = strcat('sim_',strtype,'_bird_moves_n',num2str(n_new),'_snr',num2str(snr),'_nx',num2str(nx_new),'_ny',num2str(nx_new),'_tn',num2str(tn_new),'iter',num2str(ipropnaive(niter/2)),'.mat')
    load(filename)
    subplot(3,3,6 + inx)          
    imagesc(propnaive)
    set(gca,'XTick',[], 'YTick', [])
    %xlabel(strcat('median MSE=',num2str(round(propnaiveMSE(ipropnaive(niter/2)),3))))
    
    if inx==1
        ylabel('naive')
    end

end

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'position',[pos(1),pos(2),pos(3),pos(4)*1.2])
pos = get(h,'Position');
set(h, 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(h,'PaperPositionMode','Auto')
print(h,'effect_nx','-dpdf','-r0','-bestfit')