% ProxMap for isotropic Total Variation Minimization
% Copyright: Hatef Monajemi
% 
% Date Started: 2013 July 28
% Last Modified: 2013 July 28
function [x,inform] = TVproxMap_iso(z, coef, D, dim, inform)




switch dim
 case 1
   error('Not yet implemented');
 case 2
	mu = 1/(8+.1);   % mu < 2/8 is safe for 2D. Fixed point iterations converges
 	[x , inform] = TV(z, coef, D, mu,dim, inform);
 case 3
   	mu = 1/(12+.1);   % mu < 2/12 is safe for 3D. Fixed point iterations converges
   	[x , inform] = TV(z, coef, D, mu, inform);
     
otherwise
	warning('Unexpected option');
end     
	
	
	




% Functions
function [x, inform]= TV(z, coef, D, mu ,dim, inform)

% Define an operator  H
Dz = D*z;
DDT = D*D';
%Hoperator = operator(@(v, mode) Hfunc(v, Dz, mu, coef, DDT, dim), size(D,1), size(D,1));
Hoperator_Chambolle = operator(@(v, mode) Hfunc_CBT(v, Dz, mu, coef, DDT,dim), size(D,1), size(D,1));

% Parameters for FixedPointItr
FP_crit = 10^-3;
FP_maxItr = 2000;

FP_kappa = 0;
if(~isempty(inform))
v0 =  inform.v0;
else
v0 = zeros(size(D,1),1);    
end

%%% STANDARD FPI
%[v_coef, k_FP] = FixedPointItr(Hoperator, v0, FP_kappa, FP_crit, FP_maxItr);
%%% NEWTON
%F = speye(length(v0)) - mu * DDT;
%[v_coef, k_FP] = FixedPointItr_Newton(F,Dz, coef,Hoperator, v0, FP_kappa, FP_crit, FP_maxItr);
%%% ANDERSON MIXING 
%m             = 1;
%[v_coef,k_FP] = FixedPointItr_AM(Hoperator, v0, m, FP_crit, FP_maxItr);
% [v_coef, k_FP] = FixedPointItr_BT(Hoperator, v0, FP_crit, FP_maxItr);  % anisotropic for now
% 
% x            = z - mu * D' * v_coef;


%%%%%%%%%%%%%%%%  Chambolle Projection Method %%%%%%%%%%%%%%%%%%%%%
%[v_coef, k_FP] = FixedPointItr(Hoperator_Chambolle, v0, FP_kappa, FP_crit, FP_maxItr);

[v_coef, k_FP] = FixedPointItr_BT(Hoperator_Chambolle, v0, FP_crit, FP_maxItr);  % anisotropic for now
% m              = 1;
%[v_coef, k_FP] = FixedPointItr_AM(Hoperator_Chambolle, v0, m,FP_crit, FP_maxItr) ; % anisotropic for now

%%%%%%%%%%%%%%%%  Chambolle-PDCO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[v_coef,K_FP] = projection_pdco(D, z , coef, v0, FP_crit, FP_maxItr,0);
x              = z - coef * D' * v_coef;
%fprintf('K_FP = %d\n',k_FP);


inform.v0      = v_coef;
end

             

end





function out = Hfunc(v, Dz, mu, coef, DDT, dim) 
y  = Dz + v - mu * (DDT * v) ;


m = floor(length(y)/dim);
Y = reshape(y, m, dim);

	a   = zeros(1, m);
	for i = 1:m
	a(i) = norm(Y(i,:),2);
	end


abs_y = repmat(a,dim,1);
out = y - (y./abs_y) .* proxMap(abs_y, coef/mu , 'SOFT');
end

                        
function out = Hfunc_CBT(v, Dz, mu, coef, DDT, dim)
y  = (mu/coef) * Dz + v - mu * (DDT* v) ;

m = floor(length(y)/dim);
Y = reshape(y, m, dim);

	a   = zeros(m,1);
	for i = 1:m
    a(i) = max(1,norm(Y(i,:),2));
    end

A = repmat(a,dim,1);    
out = y ./ A;    % modified version
end
           
