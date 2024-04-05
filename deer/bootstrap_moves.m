function [Xboot, yboot] = bootstrap_moves(X, y, n_replicas, n_movements, balance)

%randsample --> requires Statistics and Machine Learning Toolbox
    if balance == 1
        n_replicas_balance = [round(n_replicas/(1-mean(y))), round(n_replicas/(mean(y)))];
    else
        n_replicas_balance = [n_replicas, n_replicas];
    end
    
    Xboot = [];
    yboot = [];
    
    for i = 1:size(X,1)
        i
        x = X(i,:);
        yi = y(i);
        x_movement = [];

        for j = 1:length(x)
            x_movement = [x_movement, repmat(j,1,x(j))];
        end

        xboot = zeros(n_replicas_balance(yi+1),size(X,2));
        for nr = 1:n_replicas_balance(yi+1)
             x_movement_boot = randsample(x_movement,n_movements);
             
             for xk = x_movement_boot
                xboot(nr,xk) = xboot(nr,xk) + 1;
             end
        end

        yiboot = repmat(yi, n_replicas_balance(yi+1), 1);

        Xboot = [Xboot; xboot];
        yboot = [yboot; yiboot];

    end

