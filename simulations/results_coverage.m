n = 500;
snr = 1;
tn = 2880;
nx = 30;
ny = nx; 
maxiter = 40;

for type = 1:3
    h = figure

    if type == 1
        strtype = 'lake'
    elseif type ==2
        strtype = 'river'
    elseif type ==3
        strtype = 'sides'
    end

    prop=zeros(nx,ny);
    nx_new = nx;
    ny_new = ny; 

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
    
                if (j<i-5*ny_new/7)
                    prop(i,j)=1;
                end
            end
        end
    
    end

    msehat = zeros(1,maxiter);
    coveragehat = zeros(1,maxiter);
    widthhat = zeros(1,maxiter);
    msehat0 = zeros(1,maxiter);
    coveragehat0 = zeros(1,maxiter);
    widthhat0 = zeros(1,maxiter);
    msehat1 = zeros(1,maxiter);
    coveragehat1 = zeros(1,maxiter);
    widthhat1 = zeros(1,maxiter);

    msesmooth = zeros(1,maxiter);
    coveragesmooth = zeros(1,maxiter);
    widthsmooth = zeros(1,maxiter);
    msesmooth0 = zeros(1,maxiter);
    coveragesmooth0 = zeros(1,maxiter);
    widthsmooth0 = zeros(1,maxiter);
    msesmooth1 = zeros(1,maxiter);
    coveragesmooth1 = zeros(1,maxiter);
    widthsmooth1 = zeros(1,maxiter);

    msenaive = zeros(1,maxiter);
    coveragenaive = zeros(1,maxiter);
    widthnaive = zeros(1,maxiter);
    msenaive0 = zeros(1,maxiter);
    coveragenaive0 = zeros(1,maxiter);
    widthnaive0 = zeros(1,maxiter);
    msenaive1 = zeros(1,maxiter);
    coveragenaive1 = zeros(1,maxiter);
    widthnaive1 = zeros(1,maxiter);

    mseglm = zeros(1,maxiter);
    coverageglm = zeros(1,maxiter);
    widthglm = zeros(1,maxiter);
    mseglm0 = zeros(1,maxiter);
    coverageglm0 = zeros(1,maxiter);
    widthglm0 = zeros(1,maxiter);
    mseglm1 = zeros(1,maxiter);
    coverageglm1 = zeros(1,maxiter);
    widthglm1 = zeros(1,maxiter);

    msespline = zeros(1,maxiter);
    coveragespline = zeros(1,maxiter);
    widthspline = zeros(1,maxiter);
    msespline0 = zeros(1,maxiter);
    coveragespline0 = zeros(1,maxiter);
    widthspline0 = zeros(1,maxiter);
    msespline1 = zeros(1,maxiter);
    coveragespline1 = zeros(1,maxiter);
    widthspline1 = zeros(1,maxiter);

    iprop1 = prop==1;
    iprop0 = prop==0;

    %prop = (prop - mean(prop(:)))/std(prop(:));

    for iter = 1:maxiter
        n_movements = 720;
        n_replicas = 5000;
        nboot = 48;
        filename = strcat('BOOT_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(iter),'nboot',num2str(nboot),'tn_new', num2str(n_movements), 'n_new', num2str(n_replicas),'.mat');
        splines_filename = strcat('spline_BOOT_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(iter),'nboot',num2str(nboot),'tn_new', num2str(n_movements), 'n_new', num2str(n_replicas),'.mat');
        %filename = strcat('BOOT_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(iter),'.mat');
        load(filename)
        load(splines_filename)

        size_prop = [nx, ny];

        all_prop = all_prophat_boot;
        all_prop_vect = all_prophat_boot;

        [lower_prop, median_prop, upper_prop, lower_prop_vect, median_prop_vect, upper_prop_vect, min_val, max_val] = process_bootstrap_moves(all_prop, all_prop_vect, size_prop, prop);
        
        msehat(iter) = mean((median_prop(:)-prop(:)).^2);
        coveragehat(iter) = mean((lower_prop(:)<=prop(:)) & (upper_prop(:)>=prop(:)));
        widthhat(iter) = mean(upper_prop(:) - lower_prop(:));
        
        all_prop = all_propsmooth_boot;
        all_prop_vect = all_propsmooth_boot;

        [lower_prop, median_prop, upper_prop, lower_prop_vect, median_prop_vect, upper_prop_vect, min_val, max_val] = process_bootstrap_moves(all_prop, all_prop_vect, size_prop, prop);
        msesmooth(iter) = mean((median_prop(:)-prop(:)).^2);
        coveragesmooth(iter) = mean((lower_prop(:)<=prop(:)) & (upper_prop(:)>=prop(:)));
        coveragesmooth0(iter) = mean((lower_prop(iprop0)<=prop(iprop0)) & (upper_prop(iprop0)>=prop(iprop0)));
        coveragesmooth1(iter) = mean((lower_prop(iprop1)<=prop(iprop1)) & (upper_prop(iprop1)>=prop(iprop1)));
        widthsmooth(iter) = mean(upper_prop(:) - lower_prop(:));

        all_prop = all_propnaive_boot;
        all_prop_vect = all_propnaive_boot;

        [lower_prop, median_prop, upper_prop, lower_prop_vect, median_prop_vect, upper_prop_vect, min_val, max_val] = process_bootstrap_moves(all_prop, all_prop_vect, size_prop, prop);
        msenaive(iter) = mean((median_prop(:)-prop(:)).^2);
        coveragenaive(iter) = mean((lower_prop(:)<=prop(:)) & (upper_prop(:)>=prop(:)));
        widthnaive(iter) = mean(upper_prop(:) - lower_prop(:));

        all_prop = all_propglm_boot;
        all_prop_vect = all_propglm_boot;

        [lower_prop, median_prop, upper_prop, lower_prop_vect, median_prop_vect, upper_prop_vect, min_val, max_val] = process_bootstrap_moves(all_prop, all_prop_vect, size_prop, prop);
        mseglm(iter) = mean((median_prop(:)-prop(:)).^2);
        coverageglm(iter) = mean((lower_prop(:)<=prop(:)) & (upper_prop(:)>=prop(:)));
        widthglm(iter) = mean(upper_prop(:) - lower_prop(:));
        
        all_prop = all_propspline_boot;
        all_prop_vect = all_propspline_boot;

        [lower_prop, median_prop, upper_prop, lower_prop_vect, median_prop_vect, upper_prop_vect, min_val, max_val] = process_bootstrap_moves(all_prop, all_prop_vect, size_prop, prop);
        msespline(iter) = mean((median_prop(:)-prop(:)).^2);
        coveragespline(iter) = mean((lower_prop(:)<=prop(:)) & (upper_prop(:)>=prop(:)));
        widthspline(iter) = mean(upper_prop(:) - lower_prop(:));
        
    end
    
    mses = [mean(msehat) mean(msenaive) mean(mseglm) mean(msesmooth) mean(msespline)];
    [a, imse] = min(mses);

    coverages = [mean(coveragehat) mean(coveragenaive) mean(coverageglm) mean(coveragesmooth) mean(coveragespline)];
    [a, icoverage] = max(coverages);

    widths = [mean(widthhat) mean(widthnaive) mean(widthglm) mean(widthsmooth) mean(widthspline)];
    [a, iwidth] = min(widths);

    if type == 1
        qut_str = 'TV&';
        naive_str = 'empirical&';
        glm_str = 'GLM&';
        smooth_str = 'GPR&';
        spline_str = 'splines&';
    else
        qut_str = strcat(qut_str,'&&');
        glm_str = strcat(glm_str,'&&');
        naive_str = strcat(naive_str,'&&');
        smooth_str = strcat(smooth_str,'&&');
        spline_str = strcat(spline_str,'&&');
    end
    
    if imse == 1
        qut_str = strcat(qut_str,'{', num2str(round(mean(msehat),3)),'}');
    else
        qut_str = strcat(qut_str, num2str(round(mean(msehat),3)));
    end

    if imse == 2
        naive_str = strcat(naive_str,'{', num2str(round(mean(msenaive),3)),'}');
    else
        naive_str = strcat(naive_str, num2str(round(mean(msenaive),3)));
    end

    if imse == 3
        glm_str = strcat(glm_str,'{', num2str(round(mean(mseglm),2)),'}');
    else
        glm_str = strcat(glm_str, num2str(round(mean(mseglm),2)));
    end
    
    if imse == 4
        smooth_str = strcat(smooth_str,'{', num2str(round(mean(msesmooth),3)),'}');
    else
        smooth_str = strcat(smooth_str, num2str(round(mean(msesmooth),3)));
    end

    if imse == 5
        spline_str = strcat(spline_str,'{', num2str(round(mean(msespline),3)),'}');
    else
        spline_str = strcat(spline_str, num2str(round(mean(msespline),3)));
    end

    qut_str = strcat(qut_str,'&');
    glm_str = strcat(glm_str,'&');
    naive_str = strcat(naive_str,'&');
    smooth_str = strcat(smooth_str,'&');
    spline_str = strcat(spline_str,'&');

    if icoverage == 1
        qut_str = strcat(qut_str,'{', num2str(round(mean(coveragehat),3)),'}');
    else
        qut_str = strcat(qut_str, num2str(round(mean(coveragehat),3)));
    end

    if icoverage == 2
        naive_str = strcat(naive_str,'{', num2str(round(mean(coveragenaive),3)),'}');
    else
        naive_str = strcat(naive_str, num2str(round(mean(coveragenaive),3)));
    end

    if icoverage == 3
        glm_str = strcat(glm_str,'{', num2str(round(mean(coverageglm),3)),'}');
    else
        glm_str = strcat(glm_str, num2str(round(mean(coverageglm),3)));
    end
    
    if icoverage == 4
        smooth_str = strcat(smooth_str,'{', num2str(round(mean(coveragesmooth),3)),'}');
    else
        smooth_str = strcat(smooth_str, num2str(round(mean(coveragesmooth),3)));
    end

    if icoverage == 5
        spline_str = strcat(spline_str,'{', num2str(round(mean(coveragespline),3)),'}');
    else
        spline_str = strcat(spline_str, num2str(round(mean(coveragespline),3)));
    end
    
    qut_str = strcat(qut_str,'&');
    glm_str = strcat(glm_str,'&');
    naive_str = strcat(naive_str,'&');
    smooth_str = strcat(smooth_str,'&');
    spline_str = strcat(spline_str,'&');


    if iwidth == 1
        qut_str = strcat(qut_str,'{', num2str(round(mean(widthhat),3)),'}');
    else
        qut_str = strcat(qut_str, num2str(round(mean(widthhat),2)));
    end

    if iwidth == 2
        naive_str = strcat(naive_str,'{', num2str(round(mean(widthnaive),3)),'}');
    else
        naive_str = strcat(naive_str, num2str(round(mean(widthnaive),3)));
    end

    if iwidth == 3
        glm_str = strcat(glm_str,'{', num2str(round(mean(widthglm),3)),'}');
    else
        glm_str = strcat(glm_str, num2str(round(mean(widthglm),3)));
    end
    
    if iwidth == 4
        smooth_str = strcat(smooth_str,'{', num2str(round(mean(widthsmooth),3)),'}');
    else
        smooth_str = strcat(smooth_str, num2str(round(mean(widthsmooth),3)));
    end

    if iwidth == 5
        spline_str = strcat(spline_str,'{', num2str(round(mean(widthspline),3)),'}');
    else
        spline_str = strcat(spline_str, num2str(round(mean(widthspline),3)));
    end

    [a, iprophat] = sort(msehat);
    [mean(msehat), mean(coveragehat), mean(widthhat)]
    ihat = iprophat(round(length(iprophat)/2));

    [a, ipropsmooth] = sort(msehat);
    [mean(msesmooth), mean(coveragesmooth), mean(widthsmooth)]
    ismooth = ipropsmooth(round(length(ipropsmooth)/2));

    [a, ipropnaive] = sort(msenaive);
    [mean(msenaive), mean(coveragenaive), mean(widthnaive)]
    inaive = ipropnaive(round(length(ipropnaive)/2));

    [a, ipropglm] = sort(mseglm);
    [mean(mseglm), mean(coverageglm), mean(widthglm)]
    iglm = ipropglm(round(length(ipropglm)/2));

    [a, ipropspline] = sort(msespline);
    [mean(msespline), mean(coveragespline), mean(widthspline)]
    ispline = ipropspline(round(length(ipropspline)/2));
    
    
    itypes = [ihat, ispline, ismooth, inaive];
    for itype = 1:4

        filename = strcat('BOOT_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(itypes(itype)),'nboot',num2str(nboot),'tn_new', num2str(n_movements), 'n_new', num2str(n_replicas),'.mat');
        load(filename)

        spline_filename = strcat('spline_BOOT_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(itypes(itype)),'nboot',num2str(nboot),'tn_new', num2str(n_movements), 'n_new', num2str(n_replicas),'.mat');
        load(spline_filename)

        if itype == 1
            all_prop = all_prophat_boot;
            all_prop_vect = all_prophat_boot;
            ylab = 'TV';
        elseif itype == 5
            all_prop = all_propglm_boot;
            all_prop_vect = all_propglm_boot;
            ylab = 'GLM';
        elseif itype == 4
            all_prop = all_propnaive_boot;
            all_prop_vect = all_propnaive_boot;
            ylab = 'empirical';
        elseif itype == 3
            all_prop = all_propsmooth_boot;
            all_prop_vect = all_propsmooth_boot;
            ylab = 'GPR';
        elseif itype == 2
            all_prop = all_propspline_boot;
            all_prop_vect = all_propspline_boot;
            ylab = 'splines';
        end
        [lower_prop, median_prop, upper_prop, lower_prop_vect, median_prop_vect, upper_prop_vect, min_val, max_val] = process_bootstrap_moves(all_prop, all_prop_vect, size_prop, prop);


        subplot(4,3,3*(itype-1)+1)
        imagesc(lower_prop)
        caxis([0 1]); % Ensure consistent color limits
        set(gca,'XTick',[], 'YTick', [])
        if itype == 1
            title('lower bound')
        end
        ylabel(ylab)
    
        subplot(4,3,3*(itype-1)+2)
        imagesc(median_prop)
        caxis([0 1]); % Ensure consistent color limits
        set(gca,'XTick',[], 'YTick', [])
        if itype == 1
            title('estimate')
        end

        subplot(4,3,3*(itype-1)+3)
        imagesc(upper_prop)
        caxis([0 1]); % Ensure consistent color limits
        set(gca,'XTick',[], 'YTick', [])
        if itype == 1
            title('upper bound')
        end

    end

    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'position',[pos(1),pos(2),pos(3)*1,pos(4)*1.7])
    pos = get(h,'Position');
    set(h, 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    set(h,'PaperPositionMode','Auto')
    print(h,strcat('coverage_',strtype),'-dpdf','-r0','-bestfit')

end    

qut_str = strcat(qut_str,'\\');
naive_str = strcat(naive_str,'\\');
smooth_str = strcat(smooth_str,'\\');
glm_str = strcat(glm_str,'\\');
spline_str = strcat(spline_str,'\\');

tex_str = sprintf('%s\n',qut_str, spline_str, smooth_str, naive_str, '\hline')