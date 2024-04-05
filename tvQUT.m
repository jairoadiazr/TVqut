function sol = tvQUT(X,y,D_TV,MC,options)

    %%%%CALCULATE LAMBDA QUT
    [lambda_qut, lambda_max] = lambdaQUT(X, y, D_TV, MC);

    if lambda_qut>=lambda_max
        sol.x = zeros(size(X,2),1);
    else
		%SOLVE OPTIMIZATION
        sol = MFISTAlosses_TV(X, y, 'entropy', lambda_qut, D_TV, 2, [], options);
    end