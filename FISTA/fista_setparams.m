% Default parameters for FISTA, MFISTA, MFISTA_TV
%
% copyright 2013 Hatef Monajemi (monajemi@stanford.edu)
function defopts = fista_setparams(opts)

%-------------------------------------------
%   BEGIN DEFAULT VALUES, ALL LOWER CASE
%-------------------------------------------
% FISTA General params
defopts.fistatol		= 1e-4		; % stopping criteria for Fista
defopts.fistaitr		= 2000		; % allowable number of fista itrs
defopts.fistalog		= 1		    ; % print log
defopts.fistalogperiod  	= 10    ; % print once every 10 line 
defopts.fistax0			= []		; % initial value for FISTA
defopts.fistablocksize 	= 1			; % the block sizes for thresholding (for quaternion set this to 4); 

% Convergence studies
defopts.truesol         = []        ;
defopts.truelambda        = []      ;
defopts.tau             = 1         ; % lambda_auto (t) = tau * sigma_hat(t)


% Fixed-Point Iterations params (for TV)
defopts.fptol 		= 1e-5		; % stopping criteria for fixed-point
defopts.fpitr		= 2000		; % allowable number of fp iterations
defopts.tviso         = false     ; % if true it solves the isotropic TV
%----------------------------------------
%        END DEFAULT OPTS
%----------------------------------------

% return if no user opts
if nargin == 0 || isempty(opts) 
	return
end

% list of valid field names
vfields = fieldnames( defopts );

% Grab valid fields from user's opts
for i = 1:length(vfields)
    field = vfields{i};
    if isfield( opts, field );
        defopts.(field) = opts.(field);
    end
end



