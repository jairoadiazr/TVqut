% kappa \in [0, 1)
% H here is a proxMap operator for TV such as
% that of Chambolle, such that
% Beck and Taboule (BT) suggested
% a FISTA-type algorithm for the FPI
% and called it Fast Gradient Projection
             
function [v,k] = FixedPointItr_BT(H, v0, crit, maxItr)
	
	tol = Inf;
	k = 0;
    	s0  			  = 1;
    	y   			  = v0;
	while (tol > crit && k <= maxItr)
   	k   = k+1;
    
   	v   = (H * y)   ;
    	%=============== UPDATES =================
    	% stepsize
    	s 	= (1+ sqrt(1+ 4*s0^2))/2;
             
    	% Psudo point
    	y 	= v + ((s0-1)/s) .* (v - v0) ;
             

   	 tol = norm(v-v0, 2)/norm(v);
	% update v0 for the next itr
	v0  = v;
   	s0  = s;
	end

end

