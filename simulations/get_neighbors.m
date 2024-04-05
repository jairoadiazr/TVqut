function neighbors = get_neighbors(x, y, nx, ny)

    if x==1

        if y==1 %(1,1)
            neighbors = [2, 1; 1, 2];
        elseif y==ny %(1,ny)
            neighbors = [2, ny; 1, ny-1];
        else %(1,y)
            neighbors = [1, y-1; 1, y+1; 2, y];
        end

    elseif x==nx
        
        if y==1 %(nx,1)
            neighbors = [nx,2; nx-1, 1];
        elseif y==ny %(nx,ny)
            neighbors = [nx,ny-1; nx-1, ny];
        else %(nx,y)
            neighbors = [nx,y-1; nx-1, y; nx, y+1];
        end

    else

        if y==1 %(x,1)
            neighbors = [x, 2; x-1, 1; x+1, 1];
        elseif y==ny %(x,ny)
            neighbors = [x, ny-1; x-1, ny; x+1, ny];
        else %(x,y)
            neighbors = [x,y+1; x,y-1; x+1, y; x-1, y];
        end

    end

    neighbors = [neighbors; x,y]; %also there is the chance of not moving
end