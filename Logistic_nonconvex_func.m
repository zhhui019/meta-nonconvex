%solve the logistic regression model with the nonconvex penaty

% Inputs:
% X - a design matrix
% Y - a response vector
% maxiter  - maximal number of iterations allowed
% threshold - the nonconvex penaty type
% lambda - a tuning parameter of penalty



%
% Outputs:
% beta      - estimated coefficients
% beta0 - estimated intercept

function [beta,beta0] = Logistic_nonconvex_func(X,Y,maxiter,threshold,lambda)
                                     
    [row,col]=size(X);   
    beta_int=zeros(col,1);
    temp=sum(Y)/row;  
    beta0=log(temp/(1-temp)); 
    
   
% Step 1: Initialize (u,w,z) %
    u = exp(beta0+X * beta_int)./(1 + exp(beta0+X * beta_int));
    W = diag(u .* (1 - u));
    try
        pinv_W = pinv(W);
    catch 
        W = 1.0/4 *eye(size(W));
    end 
    z = beta0+X * beta_int + pinv(W) * (Y - u);


% Step 2: Approximate Message Passing algorithm for sparse logistic with the nonconvex penalty %
      
    iter=1;
    inter_max=1;
    beta=beta_int;
    beta_old=ones(col,1); 
    lambda_noise = 0; 
    amp0 = true;
   
 
while iter<=maxiter && norm(beta_old - beta) > (1E-4)
 %   t_start = tic;    
    beta_old=beta;
    Y_ = sqrtm(W)*z;
    X_ = sqrtm(W)*X;
    Y_ = Y_-mean(Y_);
    X_ = X_-mean(X_);
    beta0=sum(W *( z - X * beta))/sum(sum(W)); 
    beta = amp(X_/norm(X_), Y_/norm(X_), lambda_noise,amp0, inter_max,threshold,lambda); 
    X_beta = beta0 + X * beta;
    u = exp(X_beta)./(1 + exp(X_beta));
   
    W = diag(u .* (1 - u)); 
    try
        pinv_W = pinv(W);
    catch 
        W = 1.0/4 *eye(size(W));
    end
 
    z = X_beta + pinv(W) * (Y - u);

    iter=iter+1;
    
end
