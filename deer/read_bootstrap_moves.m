function [all_prop, all_prop_vect, size_prop] = read_bootstrap_moves(type, maxiter, ehdv)
 
    all_prop = [];
    all_prop_vect = [];
    
    for iter=1:maxiter
        filename = strcat('deer_EHDV-', num2str(ehdv), '_weight_1000_year_ALL_interpolation_15__PLUS_iter', num2str(iter),'.mat');
        load(filename)
        
        if type == 'h'
            prop_vect = prop_vect_hat(:);
            prop = prop_hat(:);
        elseif type == 's'
            prop_vect = prop_vect_smooth(:);
            prop = prop_smooth(:);
        elseif type == 'n'
            prop_vect = prop_vect_naive(:);
            prop = prop_naive(:);
        elseif type == 'g'
            prop_vect = prop_vect_glm(:);
            prop = prop_glm(:);
        elseif type == 'p'
            load(strcat('spline_EHDV', num2str(ehdv),'.mat'))
            prop = propspline(:,iter);
            prop = prop(:);
            proptemp = NaN(size(prop));
            proptemp(~isnan(prop_hat)) = prop(~isnan(prop_hat));
            prop = proptemp;

            prop_vect = propsplinevect(:,iter);
            prop_vect = prop_vect(:);
            %proptemp = NaN(size(prop_vect));
            %proptemp(~isnan(prophat)) = prop_vect(~isnan(prophat));
            %prop_vect = proptemp;
            
        end

        all_prop_vect = [all_prop_vect prop_vect(:)];
        all_prop = [all_prop prop(:)];

    end

    size_prop = size(prop_hat);
end