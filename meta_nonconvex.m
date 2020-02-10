% title - Solve the meta-analysis method based on nonconvex regularization problem with a single tuning parameter
% description - Jointly fit a generalized linear model with a penalty over multiple datasets. It enables heterogeneous
%               variable selections in different datasets. Fits logistic regression model.
% details - The function minimizes \eqn{-logLik + lambda * p(beta)}, where \eqn{-logLik} is the negative of the total
%           log-Likelihood from all datasets, \eqn{lambda} is a single tuning parameter and \eqn{p(beta)} is a specific penalty
%           function enabling heterogeneous selections of variables in different datasets. For more details of the penalty
%           function, see the reference below.

% % Inputs:
% X - a concatenated design matrix, of dimension \eqn{N * p}, where N is the total sample size over
%    multiple datasets and p is the total number of variables.
% Y - a concatenated response vector from all datasets
% obs - a vector of sample sizes of multiple datasets
% lambda - a tuning parameter of penalty
% maxit  - maximal number of iterations allowed
% threshold - the nonconvex penaty type
% tol - tolerance level of convergence

% Outputs:
% mycoef - estimated coefficients in each dataset 
% intercept - estimated intercept in each dataset
% iteration - number of iterations converge, TRUE if convergence is achieved
% diff -last step difference
% @references Hui Zhang, Shou-Jiang Li, Hai Zhang, Zi-Yi Yang, Yan-Qiong Ren, Liang-Yong Xia and Yong Liang,
%     Meta-Analysis Based on Nonconvex Regularization.

% Author:   zhanghui
% email:    zhanghui.nwu@foxmail.com
% Website:  
% Last revision: 01-March-2019


function [mycoef,intercept,iteration,converge,diff] = meta_nonconvex(X,Y, obs,lambda, maxit,threshold,tol)

%starting and ending index of each dataset
   start_idx =cumsum(obs) + 1 - obs; % starting index of each dataset
   end_idx  = cumsum(obs) ; %ending index of each dataset
   M = length(obs);              %number of datasets
   p = size(X,2);                 %number of covariates
   N = sum(obs);                %total number of obserations
   gamma = zeros(p,maxit+1);     %iterations of gamma
   gamma(:,1) = ones(p,1);
   theta = zeros(M*p,maxit);      %iterations of theta
   X_tha = zeros(N,p);          % colMultiply(X.all, theta)
   beta_hat = zeros(p,maxit,M);       %iterations of beta.hat
   beta0_hat = zeros(maxit,M);
   mycoef = zeros(p,M);        %final estimate of coefficients
   intercept = zeros(1,M);     %final estimate of intercept
   itr = 1;
   m_diff = zeros(1,M);         %marginal error
   
   while itr <= maxit
       for m = 1:M
          Y_mm = Y(start_idx(m):end_idx(m));
          Y_m = Y_mm;
          X_mm = X(start_idx(m):end_idx(m),:).*gamma(:,itr)';
          X_m = X_mm;    
       
          theta_fit_beta = Logistic_nonconvex_func(X_m,Y_m,maxit,threshold,lambda);
          

          theta(((m - 1) * p + 1):(m * p), itr) = theta_fit_beta;
           %adjust X.all by colMultiply(X.all, theta) for further usage
           X_tha(start_idx(m):end_idx(m),:)=X(start_idx(m):end_idx(m),:).*theta_fit_beta';
           beta_hat(:,itr,m) = theta(((m - 1) * p + 1):(m * p), itr).*gamma(:,itr);
           beta0_hat(itr,m) = Newton_Raphson_agorithm(X(start_idx(m):end_idx(m),:),Y_m,beta_hat(:,itr,m));
           %calculate iteration difference
           if itr == 1
               m_diff(m) = max(abs(beta_hat(:,itr,m)));
           else
               m_diff(m) = max(abs(beta_hat(:,itr,m)-beta_hat(:,itr-1,m)));
           end
       end
           if max(m_diff) < tol     %break iterations if diff < tol
               break
           end
               
            itr = itr+1;            
            gamm_fit_beta = Logistic_nonconvex_func(X_tha,Y,maxit,threshold,lambda);
            gamma(:,itr) =  gamm_fit_beta;
   end
   
       
   %determine if convergence is achieved
   if  itr == 1
       iteration = itr;
       converge = false;
   elseif itr > maxit
       iteration = itr-1;
       converge = false;
   else 
       iteration = itr;
       converge = true;
   end
  
   for m=1:M
       mycoef(:,m) = beta_hat(:,iteration,m);
       intercept(m) = beta0_hat(iteration,m);
   end
   
   diff = max(m_diff);
       
            
       
           