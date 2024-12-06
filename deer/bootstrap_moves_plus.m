function [Xboot, yboot] = bootstrap_moves_plus(X, y, n_replicas, n_movements)
    
    flag = true;

    while flag == true
        selected_i = randsample(1:length(y), n_replicas, true);
        if std(y(selected_i))==0
            flag = true;
        else
            flag = false;
        end
    end

    Xboot = [];
    yboot = [];
    for i = selected_i

        x = X(i,:);
        yi = y(i);
        x_movement = [];
        
        for j = 1:length(x)
            x_movement = [x_movement, repmat(j,1,x(j))];
        end
        
        xboot = zeros(1,size(X,2));
        x_movement_boot = randsample(x_movement,n_movements);
             
        for xk = x_movement_boot
            xboot(xk) = xboot(xk) + 1;
        end
        

        yiboot = yi;

        Xboot = [Xboot; xboot];
        yboot = [yboot; yiboot];

    end
end

