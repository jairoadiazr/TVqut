
nx=50; %x grid
ny=nx; %y grid
tn = 30*24*4; %time spent per bird --> 1 month, sampled every 15 minutes
n=5000;

% mu=nx*ny propensity matrix

prop=zeros(nx,ny);

if type == 1
    strtype = 'lake'

    for i = 1:nx
        for j = 1:ny
            if (i-4*ny/7)^2 + (j-4*ny/7)^2 <= (ny/5)^2
                prop(i,j)=1;
            end
        end
    end

elseif type == 2

    strtype = 'river'
%prop=zeros(nx,ny);

   for i = 1:nx
        for j = 1:ny
            if i>(2*j-3*ny/5)&&(i<2*j-ny/6)
                prop(i,j)=1;
            end
        end
   end
   %mean(prop(:))
   %imagesc(prop)

elseif type ==3

    strtype = 'sides'

    for i = 1:nx
        for j = 1:ny
            if (i-ny/3)^2 + (j-2*ny/3)^2 <= (ny/5)^2
                prop(i,j)=1;
            end

            %if j>(i-2*ny/3)&&(j<i-4*ny/7)
            if (j<i-5*ny/7)
                prop(i,j)=1;
            end
        end
    end

end    

%bird movements
M = zeros(n, tn);
cells = reshape(1:(nx*ny),nx,ny);

%'calculating movements'
parfor i = 1:n 
    i
    %make one bird to move
    cell_movement = zeros(1,tn);
    mx0 = randsample(nx, 1);
    my0 = randsample(ny, 1);
    cell_movement(1) = cells(mx0,my0);

    for j = 1:(tn-1)
        neighbors = get_neighbors(mx0, my0, nx, ny);

        weights = ones(1,size(neighbors,1));
        for k = 1:size(neighbors,1)
            if prop(neighbors(k,1), neighbors(k,2)) > 0
                weights(k) = 2;
            end
        end
        if i>(n/2)
            weights = 3-weights;
        end
        mxy = neighbors(randsample(size(neighbors,1),1,true,weights),:);
        
        cell_movement(j+1) = cells(mxy(1), mxy(2));
        mx0 = mxy(1);
        my0 = mxy(2);
    end

    M(i,:) = cell_movement;

end

X = zeros(n, nx*ny);

parfor i = 1:n 
    %make one bird to move
    movement = zeros(1,nx*ny);
    for j = 1:tn
        movement(M(i,j)) = movement(M(i,j)) + 1;
    end
    X(i,:) = movement;
end

%save big simulations
filename = strcat('test_',strtype,'_bird_moves_main.mat')
save(filename,'prop','M','X','n','nx','ny','tn','-v7.3')

