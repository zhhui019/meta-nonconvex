%solve the intercept by Newton Raphson agorithm

function beta0_hat = Newton_Raphson_agorithm(X,Y,beta)
    [row,col]=size(X);
    temp=sum(Y)/row;  
    beta0=log(temp/(1-temp)); 
    beta_int = [beta0;beta];
% Step 1: Initialize (u,w,z) %
    X = [ones(row,1),X];
    u = exp(X * beta_int)./(1 + exp(X * beta_int));
    W = diag(u .* (1 - u));
    try
        pinv_W = pinv(W);
    catch 
        W = 1.0/4 *eye(size(W));
    end
    z = X * beta_int + pinv(W) * (Y - u);
    
    iter=1;
    maxiter=100;
    beta_hat=beta_int;
    beta_old=ones(col+1,1);
          
while iter<=maxiter && norm(beta_old - beta_hat) > (1E-8)
    beta_old=beta_hat;
    beta_zero=sum(W *( z - X(:,2:end) * beta_hat(2:end)))/sum(sum(W)); 
    beta_hat(1)=beta_zero;
    u = exp(X * beta_hat)./(1 + exp(X * beta_hat));
    W = diag(u .* (1 - u));
    try
        pinv_W = pinv(W);
    catch 
        W = 1.0/4 *eye(size(W));
    end
    z = X * beta_hat + pinv(W) * (Y - u);
    iter=iter+1;
end
beta0_hat = beta_hat(1);
