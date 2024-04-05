%reshape movement parameters from main simulation



%nx_new = 10; %discretization
%ny_new = nx_new;
%n_new = 100; %number of birds
%tn_new = nx_new*ny_new; %number of samples

cells = reshape(1:(nx*ny),nx,ny);
cells_new = reshape(1:(nx_new*ny_new),nx_new,ny_new);
cells_inter = imresize(cells_new, [nx,ny], 'nearest');

M_new = cells_inter(M);

%clear M

M_new = M_new(round(linspace(1,n,n_new)), round(linspace(1,tn,tn_new)));
y_new = y(round(linspace(1,n,n_new)));
p_new = p(round(linspace(1,n,n_new)));

clear y
clear p

propvec_new = propvec(round(linspace(1,nx,nx_new)),round(linspace(1,ny,ny_new)));

clear propvec

%bird movements
X = zeros(n_new, nx_new*ny_new);

%'calculating movements'
parfor i = 1:n_new 
   %make one bird to move
   movement = zeros(1,nx_new*ny_new);
   for j = 1:tn_new
       movement(M_new(i,j)) = movement(M_new(i,j)) + 1;
   end
   X(i,:) = movement;
end

clear M_new

nx = nx_new;
ny = ny_new;
n = n_new;


%Calculate D matrix
diagonal = 0;

if diagonal==1
    D_TV = zeros(3*(nx-1)*ny+(ny-1)*nx, nx*ny);
else
    D_TV = zeros((nx-1)*ny+(ny-1)*nx, nx*ny);
end


for(i=1:1:(nx-1))
    for(j=1:1:ny)
        D_TV((j-1)*(nx-1)+i,(j-1)*(nx-1)+j-1+(i:(i+1)))=[-1 1];
    end
end
k0=(nx-1)*ny;
for(i=1:1:(ny-1))
    for(j=1:1:nx)
        D_TV(k0+(i-1)*nx+j,(i-1)*nx+[j j+nx]  )=[-1 1];
    end
end

if diagonal==1
    k0=2*(nx-1)*ny;
    for(i=1:1:(ny-1))
        for(j=1:1:(nx-1))
            D_TV(k0+(i-1)*nx+j,(i-1)*nx+[j j+nx+1]  )=[-1 1];
        end
    end
    
    k0=3*(nx-1)*ny;
    for(i=1:1:(ny-1))
        for(j=2:1:nx)
            D_TV(k0+(i-1)*nx+j,(i-1)*nx+[j j+nx-1]  )=[-1 1];
        end
    end
    
    D_TV = D_TV(sum(D_TV'~=0)~=0,:);
end

D_TV=sparse(D_TV);

y = y_new;
p = p_new;
propvec = propvec_new;


