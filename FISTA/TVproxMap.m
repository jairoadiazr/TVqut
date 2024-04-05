% ProxMap for Total Variation Minimization
% Copyright: Hatef Monajemi
% 
% Date: 2013 July 17
function [x,inform] = TVproxMap(z, coef, D, dim, inform)




switch dim
 case 1
        mu = 1/(4+.1);   % mu < 2/8 is safe for 2D. Fixed point iterations converges
 	[x , inform] = TV(z, coef, D, mu, inform);
 case 2
	mu = 1/(8+.1);   % mu < 2/8 is safe for 2D. Fixed point iterations converges
 	[x , inform] = TV(z, coef, D, mu, inform);
 case 3
   	mu = 1/(12+.1);   % mu < 2/12 is safe for 3D. Fixed point iterations converges
   	[x , inform] = TV(z, coef, D, mu, inform);
     
otherwise
	warning('Unexpected option');
end     
	
	
	




% Functions
function [x, inform]= TV(z, coef, D, mu , inform)

% Define an operator  H
Dz = D*z;
DDT = D*D';
Hoperator = operator(@(v, mode) Hfunc(v, Dz, mu, coef, DDT), size(D,1), size(D,1));
Hoperator_Chambolle = operator(@(v, mode) Hfunc_CBT(v, Dz, mu, coef, DDT), size(D,1), size(D,1));

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







function out = Hfunc(v, Dz, mu, coef, DDT) 
y  = Dz + v - mu * (DDT * v) ;
out = y - proxMap(y, coef/mu , 'SOFT');
end

                        
function out = Hfunc_CBT(v, Dz, mu, coef, DDT)
y  = (mu/coef) * Dz + v - mu * (DDT* v) ;
out =  y ./ max(1, abs(y));      % Modified
%out  = y ./ ( 1 + abs(y-v) )  ; % Original 
end
                        

end
