NXs = [10];
Ns = 50:50:400;
snr = 1;
niter = 400;
tn = 30;
strtype = 'lake';


P = NaN(4, length(Ns));
P0 = NaN(4, length(Ns));

for inx = 1:length(NXs)
    nx = NXs(inx);
    ny = nx;

    for in = 1:length(Ns)
        n = Ns(in);
        %[nx, n]

        test = zeros(1,niter);
        test_lasso = zeros(1,niter);
        test_llratio = zeros(1,niter);
        test_llratio_exact = zeros(1,niter);

        for iter = 1:niter
            clear lambda_max
            clear lambda_lasso_max
            clear lambda_qut
            clear lambda_lasso_qut

            filename = strcat('test_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(iter),'.mat');
            load(filename)
            filename
            [lambda_qut, lambda_max];
            test(iter) = lambda_qut<lambda_max;  
            test_lasso(iter) = lambda_lasso_qut<lambda_lasso_max;  

            if n>(nx*ny)
                test_llratio(iter) = llratio;
                test_llratio_exact(iter) = llratio_exact;
            end
        end
        P(1,in) = mean(test);
        P(2,in) = mean(test_lasso);

        if n>(nx*ny)
            P(3,in) = mean(test_llratio);
            P(4,in) = mean(test_llratio_exact);
        end

    end
end


snr = 0;
for inx = 1:length(NXs)
    nx = NXs(inx);
    ny = nx;

    for in = 1:length(Ns)
        n = Ns(in);
        %[nx, n]

        test = zeros(1,niter);
        test_lasso = zeros(1,niter);
        test_llratio = zeros(1,niter);
        test_llratio_exact = zeros(1,niter);

        for iter = 1:niter
            filename = strcat('test_sim_',strtype,'_bird_moves_n',num2str(n),'_snr',num2str(snr),'_nx',num2str(nx),'_ny',num2str(ny),'_tn',num2str(tn),'iter',num2str(iter),'.mat')
            load(filename)
            test(iter) = lambda_qut<lambda_max;    
            test_lasso(iter) = lambda_lasso_qut<lambda_lasso_max; 

            if n>(nx*ny)
                test_llratio(iter) = llratio;
                test_llratio_exact(iter) = llratio_exact;
            end
        end
        P0(1,in) = mean(test);
        P0(2,in) = mean(test_lasso);

        if n>(nx*ny)
            P0(3,in) = mean(test_llratio);
            P0(4,in) = mean(test_llratio_exact);
        end

    end
end

h= figure

subplot(1,2,1)
plot(Ns,P0)
yline(0.05, '--')
xlabel('n')
title('level')
xlim([50,400])

legend('TV-test', 'LASSO-test', 'likelihood ratio test', 'exact likelihood ratio test','location','west')

subplot(1,2,2)
plot(Ns,P)
xlabel('n')
title('power')
xlim([50,400])


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'position',[pos(1)*1,pos(2)*1,pos(3)*0.8,pos(4)*0.8])
pos = get(h,'Position');
set(h, 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(h,'PaperPositionMode','Auto')
print(h,'level_power','-dpdf','-r0','-bestfit')