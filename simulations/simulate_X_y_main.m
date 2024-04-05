filename = strcat(strtype,'_bird_moves_main','_iter',num2str(iter),'.mat')

load(filename)

X = zeros(n, nx*ny);

'calculating X matrix'
parfor i = 1:n 
    %make one bird to move
    movement = zeros(1,nx*ny);
    for j = 1:tn
        movement(M(i,j)) = movement(M(i,j)) + 1;
    end
    X(i,:) = movement;
end

'calculating y'

propvec = reshape(snr*prop/100, nx*ny, 1);
mu=X * propvec;
p=1./(1+exp(-mu));
y=binornd(1,p);

mean(y)

propvec = reshape(propvec,nx,ny);

%save big simulations
clear X
filename = strcat(strtype,'_bird_moves_main_snr',num2str(snr),'iter',num2str(iter),'.mat');
save(filename,'M','y','p','propvec','nx','ny','tn','n','snr','-v7.3')



