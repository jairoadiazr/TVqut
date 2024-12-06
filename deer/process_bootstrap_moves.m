function [lower_prop, median_prop, upper_prop, lower_prop_vect, median_prop_vect, upper_prop_vect, min_val, max_val] = process_bootstrap_moves(all_prop, all_prop_vect, size_prop, prop)
    
    prop_boot = mean(all_prop');
    all_prop_bc = 2*all_prop-prop_boot';
    lower_prop = reshape(quantile(all_prop_bc',0.025), size_prop);
    median_prop = reshape(mean(all_prop_bc'), size_prop);
    upper_prop = reshape(quantile(all_prop_bc',0.975), size_prop);

    prop_vect_boot = mean(all_prop_vect');
    all_prop_vect_bc = 2*all_prop_vect-prop_vect_boot';
    lower_prop_vect = quantile(all_prop_vect_bc',0.025);
    median_prop_vect = mean(all_prop_vect_bc');
    upper_prop_vect = quantile(all_prop_vect_bc',0.975);
    
    min_val = min(lower_prop(:));
    max_val = max(upper_prop(:));
    
    %mean((prop(:)-median_prop(:)).^2);
     
    if isempty(prop)==0
        %A = [median_prop(:), ones(length(median_prop(:)),1)];
        %A_lower = [lower_prop(:), ones(length(median_prop(:)),1)];
        %A_upper = [upper_prop(:), ones(length(median_prop(:)),1)];
        %Y = prop;
        %coeffs = A \ Y(:);
        %median_prop = reshape(A*coeffs, size_prop);
        %lower_prop = reshape(A_lower*coeffs, size_prop);
        %upper_prop = reshape(A_upper*coeffs, size_prop);
        %mean((prop(:)-median_prop(:)).^2);

        % minprop = min(median_prop(:));
        % median_prop = median_prop - minprop;
        % maxprop = max(median_prop(:));
        % median_prop = median_prop/maxprop;
        % 
        % upper_prop = (upper_prop-minprop)/maxprop;
        % lower_prop = (lower_prop-minprop)/maxprop;

        minprop = min(upper_prop(:));
        tmp_prop = lower_prop - minprop;
        maxprop = max(tmp_prop(:));
        
        median_prop = median_prop-min(median_prop(:));
        median_prop = median_prop/max(median_prop(:));

        upper_prop = (upper_prop-minprop)/maxprop;
        lower_prop = (lower_prop-minprop)/maxprop;
        
        %upper_prop = (upper_prop - mean(median_prop(:)))/std(median_prop(:));
        %lower_prop = (lower_prop - mean(median_prop(:)))/std(median_prop(:));
        %median_prop = (median_prop - mean(median_prop(:)))/std(median_prop(:));

        
    end
end