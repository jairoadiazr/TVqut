function [w, lambda] = SolveInfinityNormDeer_cvx(y,X,D)

m 		= size(D,1);
X0 = X*ones(size(X,2),1);
beta0 = glmfit(X0,y,'binomial','Constant','off');
u0 = beta0*ones(size(X,2),1);
p = 1./(1+exp(-X*u0));
grad_h = X'*(y-p);

cvx_begin quiet
variable w(m);
minimize norm(w,inf)
subject to
   D' * w == grad_h
cvx_end
lambda = cvx_optval;

end

