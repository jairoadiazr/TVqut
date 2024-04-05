% This routine solves the TV minimization
% problem using Monotone Fast Iterative Shrinkage
% Tresholding Algorithm (MFISTA) 
% 
% Copyright 2012, Hatef Monajemi (monajemi@stanford.edu)
% http://www.stanford.edu/~monajemi

% Date Started: 22 Feb 2012
% Last Modified:03 Feb 2014
%               12 MAY 2014 (merged the iso part with this code)

% input : 
% A : an n-by-m explicit matrix or an operator. An operator (refer to ASP code of Michael Saunders and see '@operator')
% example:
% A  = operator(@(x,mode) myfunc(x, mode),m,n);
% If mode == 1, returns  v1 =   A * x;
% if mode == 2, returns  v2 =   A'* x; 
%
% b 	:  data 
% lambda:  regularization parameter
% D     :  derivative matrix
% dim	:  dimnsion of problem 1|2|3
% Lf    :  Lf (Lipschitz constant)   ( [] : equivalent to backtracking )
% options: options for FISTA and Fixed-Point. Refer to fista_setparams.m
%
%
% \frac{|F(x_k) - mu_F(x_k)|}{mu_F(x_k)} < crit,
% where 'mu_F' is the 'average' of min{k,10} prior values of F(x)
% This termination criterion is based on NESTA paper (Becker, Bobin, Candes)
%
% output
% sol.x : argument that minimizes F(x) = f(x) + g(x) when 
% 	  anisotropic TV: f(x) = 1/2 || Ax-b ||_2^2 and g = lambda * || D x ||_1
%     isotropic TV: f(x) = 1/2 || Ax-b ||_2^2   and g = lambda * || D x ||_2,1
% sol.w : dual solution
% sol.k : number of Fista iterations to solve the problem


function sol = MFISTAlosses_TV(A, b, losstype, lambda, D, dim, Lf, options) 

% ---------------------------------
% Check number of input args
% ---------------------------------
if (nargin < 6 || nargin > 8)
     error('Wrong number of arguments.');
end

if (nargin > 5 && ~isempty(Lf) )
	back_track = false;
else
	back_track = true;
end


% Set default parameters or grab user-defined params
if (nargin < 7 || isempty(options))
     opts = fista_setparams();
else
     opts = fista_setparams(options);
end


%-----------------------------------------
%         Extract information from opts
%-----------------------------------------

%x0   	= opts.fistax0;
x0   	= options.x0;

k_max	= opts.fistaitr;
qPrint	= opts.fistalog;
qPeriod = opts.fistalogperiod;
crit    = opts.fistatol;
qIso    = opts.tviso;


% Set x0 to zero if empty
if (isempty(x0)) 
    x0 = zeros(size(A,2),1);
end


%-----------------------------------------
% define objective function F and 
% quadratic approximation of the objective 
% function F(x) at a given point Y
%-----------------------------------------
if(qIso==true)
    if(strcmp(losstype, 'l2'))
        F_Func = @(A,x,b, lambda, dim) 0.5* (norm(A*x - b, 2))^2 + lambda* BlockOneNorm(D*x,dim);
        Q_Func = @(A,x,b,lambda,y,Lp, dim) 0.5* (norm(A*y - b, 2))^2 + (x-y)' *  ( A'*(A*y-b) ) ...
                        + 0.5*Lp*(norm(x-y,2))^2  + lambda* BlockOneNorm(D*x,dim);
    end
else
    if(strcmp(losstype, 'l2'))
        F_Func = @(A,x,b, lambda,dim) 0.5* (norm(A*x - b, 2))^2 + lambda* norm(D * x,1);
        Q_Func = @(A,x,b,lambda,y,Lp, dim) 0.5* (norm(A*y - b, 2))^2 + (x-y)' *  ( A'*(A*y-b) ) ...
                              + 0.5*Lp*(norm(x-y,2))^2  + lambda* norm(D*x,1);
    end
    if(strcmp(losstype, 'sqrtl2'))
        F_Func = @(A,x,b, lambda,dim) norm(A*x - b, 2) + lambda* norm(D * x,1);
        Q_Func = @(A,x,b,lambda,y,Lp, dim) norm(A*y - b, 2) + (x-y)' *  ( A'*(A*y-b) )/norm(A*y-b,2) ...
                              + 0.5*Lp*(norm(x-y,2))^2  + lambda* norm(D*x,1);
    end
    if(strcmp(losstype, 'entropy'))
        F_Func = @(A,x,b, lambda,dim) sum(A*x .* (1-b)) + sum(log(1+exp(-A*x))) + lambda* norm(D * x,1);
        Q_Func = @(A,x,b,lambda,y,Lp, dim) sum(A*y .* (1-b)) + sum(log(1+exp(-A*y))) + (x-y)' *  ( A'*((1-b)-1./(1+exp(A*y))) ) ...
                              + 0.5*Lp*(norm(x-y,2))^2  + lambda* norm(D*x,1);
    end       
end





tol               = Inf;
k 			  	  = 0;
s0  			  = 1;
y   			  = x0;
prior_FuncVals    = zeros(1,10);
prior_FuncVals(1) = F_Func(A,x0,b, lambda,dim);
F0                = prior_FuncVals(1);

prior_FuncVals_100    = zeros(1,100);
prior_FuncVals_100(1) = F_Func(A,x0,b, lambda,dim);


%-----------------------------------------------------------------------
% 		Print log header.
%-----------------------------------------------------------------------
    if (qPrint==1) 
      fprintf('\n');
      fprintf(' %s\n',repmat('=',1,80));
      fprintf(' MFISTA_TV  (%s)\n', date);
      fprintf(' Copyright 2012, Hatef Monajemi (monajemi@stanford.edu)\n');
      fprintf(' %s\n',repmat('=',1,80));
      fprintf(' %-20s: %8i %5s'    ,'dim'  	            , dim   , '');
      fprintf(' %-20s: %8.2e\n'  	   ,'lambda'                , lambda);
      fprintf(' %-20s: %8i %5s'    ,'No. rows'          ,size(A,1)       ,'');
      fprintf(' %-20s: %8i \n'    ,'No. columns'       ,size(A,2)    );
      fprintf(' %-20s: %8.2e %5s'    ,'FISTA Optim. tol'    ,crit ,''    );
      fprintf(' %-20s: %8i \n'    ,'FISTA Max itr',k_max);
      fprintf(' %-20s: %8.2e %5s'    ,'Fixed-Point tol'       ,opts.fptol ,'');
      fprintf(' %-20s: %8i \n'    ,'Fixed-Point Max itr',opts.fpitr);
      fprintf(' %s\n',repmat('=',1,80));
      fprintf('%-20s %-20s %-20s %-20s \n', 'itr #','rel. err','func. value', 'Lf')
     end


if (back_track==false)
%===========================
% FISTA WITHOUT BACKTRACKING
%===========================
    inform = [];  %% for warm start of fixed-point iterations in prox-map.
	while (tol > crit && k<= k_max)
  
 		k 		= k+1; 
 		z 		= A*y - b;
 		tmp 	= A' * z;
 		tmp 	= y - tmp./Lf;
        %-------------
 		% Tresholding
        %-------------
        if(qIso==true)
        [x1, inform] = TVproxMap_iso(tmp, lambda/L_bar, D, dim, inform);
        else
 		[x1,inform] = TVproxMap(tmp, lambda/Lf, D, dim, inform);
        end

        %------------------------------------------
        %               UPDATES
        %------------------------------------------
        s1 		= (1+ sqrt(1+ 4*s0^2))/2; % stepsize

        %-----------------------------------------------------------------
		% Find argmin{F(x): x = x_{k-1}, x1}  (ONLY DIFFERENCE WITH FISTA)
        %-----------------------------------------------------------------
        F_k  	= F_Func(A,x1,b,lambda,dim);
        if( F_k <  F0 )
 		y 		= x1 + ((s0-1)/s1) .* (x1 - x0) ;   % Psudo point
        else
		y       = x0 + (s0/s1) .* (x1 - x0);        % Psudo point
		end

 		F_mu 	= sum(prior_FuncVals) / nnz(prior_FuncVals);
 		tol 	= abs( F_k - F_mu) / F_mu;          % tolerance
 		TOL(k) = tol;
        %-----------------------------------------------------------
 		% prior values of objective function for termination purpose
        %-----------------------------------------------------------
 		prior_FuncVals(1 + mod(k,10)) = F_k;
 
 		% update varaibles for the next iteration
		x0 	= x1;
 		s0 	= s1;  
        F0      = F_k;
 	 
        if(mod(k,qPeriod)==0 & qPrint==1)
 	    	fprintf('%-20i %-20e %-20e %-20f \n', k, tol, F_k, Lf); 
        end
	end
	
else
% FISTA WITH BACKTRACKING
    %Lf = 1;
    %Lf   = 131072;
    Lf = 524288;
    %Lf = 524288/2;
    %Lf = 262144*1024;
    'BE CAREFUL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    eta  = 2.0;
    inform = [];    
	while ( (tol > crit && k<= k_max) )

 		k = k+1;
 		if(strcmp(losstype, 'l2'))
            z 		= A  * y - b;  
            tmp1 	= A' * z;
        end
        if(strcmp(losstype, 'entropy'))
            tmp1 = A'*((1-b)-1./(1+exp(A*y)));
        end
            %-------------------------
 			%  backtracking for step k
            %-------------------------
 			ik   = 0;
 			beta = Inf;
 			while(beta > 0)
                L_bar   = (eta^ik ) * Lf;
                tmp 	= y - tmp1./L_bar;
                %-------------
                % Tresholding
                %-------------
                if(qIso==true)
                    [x1, inform] = TVproxMap_iso(tmp, lambda/L_bar, D, dim, inform);
                else
                    [x1,inform] = TVproxMap(tmp, lambda/Lf, D, dim, inform);
                end
                %-------------------------------
                % update backtracking parameters
                %-------------------------------
                beta 	= real(F_Func (A, x1, b, lambda,dim) - Q_Func(A, x1, b, lambda, y, L_bar,dim));
                ik 		= ik+1;
			end
			
                
                    
            %------------------------------------------
            %               UPDATES
            %------------------------------------------
            Lf 		=  L_bar;
 			s1 		= (1+ sqrt(1+ 4*s0^2))/2;       % stepsize
 
            %-----------------------------------------------------------------
            % Find argmin{F(x): x = x_{k-1}, x1}  (ONLY DIFFERENCE WITH FISTA)
            %-----------------------------------------------------------------
            F_k  	= F_Func(A,x1,b,lambda,dim);
            if( F_k <=  F0)
            y 		= x1 + ((s0-1)/s1) .* (x1 - x0) ; % Psudo point
            else
            y      = x0 + (s0/s1) * (x1 - x0)       ; % Psudo point
            end


 		
 		F_mu 	= sum(prior_FuncVals) / nnz(prior_FuncVals);
 		
 		tol 	= abs( F_k - F_mu) / F_mu;             % tolerance
 		TOL(k) = tol;
        %-----------------------------------------------------------
        % prior values of objective function for termination purpose
        %-----------------------------------------------------------
 		prior_FuncVals(1 + mod(k,10)) = F_k;
        
        F_mu_100 	= sum(prior_FuncVals_100) / nnz(prior_FuncVals_100);
        prior_FuncVals_100(1 + mod(k,100)) = F_k;
        if F_mu_100<F_k
            Lf = Lf*2;
        end
 		% update varaibles for the next iteration
		x0 	= x1;
 		s0 	= s1;  
        F0      = F_k;
        if(mod(k,qPeriod)==0 & qPrint==1)
 		    fprintf('%-20i %-20e %-20e %-20f \n', k, tol, F_k, Lf); 
        end
	
      end

end


sol.x  = x1;
sol.k  = k ;
sol.w  = (inform.v0)*lambda;  

if(qPrint==1)
if( k <= k_max && tol < crit)
  fprintf('FISTA: converged in %i iterations (nnz(x)= %i)\n', k, nnz(x1));
  fprintf('FISTA: Relative error: %-.2e\n', tol); 
else
  warning(['FISTA: not converged to desired tolerance, relative error achieved: ', num2str(tol),'( nnz(x)= ', num2str(nnz(x1)) ,')']);
end
end


%figure
%semilogy(TOL, 'linewidth', 3)
%set(gca, 'fontsize', 17)
%ylabel('rel. err','fontsize', 18)
%xlabel('iteration','fontsize', 18)
%title('Convergence - FISTA','fontsize', 18)
end

%--------------------------------
% Block one norm for isotropic TV
%--------------------------------
function f = BlockOneNorm(y,dim)
m = length(y)/dim;
B = reshape(y, m, dim);
a   = zeros(1, m);
for i = 1:m
        a(i) = norm(B(i,:),2);
end
f = norm(a,1);
end
 
