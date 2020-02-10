function x_pre = amp(A, y, lambda_noise, amp0, inter_max,threshold,lambda)


% Inputs:
% A         - measurement matrix
% y         - measurement vector
% x_true    - the uknown vector (here used to get mse accross iterations)
% amp0      - boolean indicator (true for noiseless case)
% inter_max - maximum number of amp iterations
% threshold - the nonconvex penaty type
% lambda    - denoiser parameter in the noisy case

%
% Outputs:
% x         - the last estimate
% tau2      - vector of effective variances accross interations
% mse       - vector of mean squared errors accross interations

% @references Hui, Zhang, Hai, et al. Approximate Message Passing Algorithm for Nonconvex Regularization[J]. IEEE Access, 2019.
%             张会, 张海. 基于AMP的L_(1/2)正则化方法[J]. 中国科学:信息科学, 2017(01):62-76.

% Author:   zhanghui
% email:    zhanghui.nwu@foxmail.com
% Website:  
% Last revision: 01-Aug-2017

[M, N] = size(A);
delta = M/N;
tau2 = zeros(inter_max, 1);

z = y;
x_hat = zeros(N,1);
tau2(1) = 1/M*(norm(z,2)^2);

p=2; %mcp parameter
a=3.7; %scad parameter
alpha = 1; %half parameter
inter_max = inter_max+1;
for t=2:inter_max
    
    % depending on if there is noise in the measuerements (amp0) the denoiser
    % parameter is calculated differenetly
    if(amp0 == true)
        deniser_parameter = sqrt(tau2(t-1))*lambda;
    else
        deniser_parameter = sqrt(lambda_noise+lambda_noise*tau2(t-1));
    end
 switch  threshold
     case 'lasso'

        x_hat = wthresh(A'*z +x_hat,'s', deniser_parameter);
        eta_prime = 1/N * nnz(wthresh(A'*z + x_hat,'s', deniser_parameter));
        z = y - A*x_hat + 1/delta * z * eta_prime;

     case 'mcp'
        
        b = A'*z + x_hat;
        for i = 1:N
            if abs(b(i)) > p*deniser_parameter            
                x_hat(i) = b(i) ;
            else
                if abs(b(i)) < p*deniser_parameter && abs(b(i)) > deniser_parameter
                x_hat(i) = (b(i) - sign(b(i))*deniser_parameter)*(p/(p-1));
                else
                    x_hat(i) = 0;
                end
            end
        end

        z = y - A*x_hat +  1/delta * z*1/N *sum(etaprimemcp(b, deniser_parameter,p));
   
      case 'half'

        b = alpha*A'*z + x_hat;
        for i = 1:N
            if abs(b(i)) >= deniser_parameter
                phi = acos(deniser_parameter/8*((abs(b(i))/3)^(-1.5)));
                x_hat(i)=4*b(i)/3*(cos(((pi)/3)-phi/3))^2;
            else
                x_hat(i) = 0;
            end
        end 

        z = y - A*x_hat +  1/delta * z*1/N *sum(etaprimehalf(b, deniser_parameter));
   
     case 'scad'
      
      b = A'*z + x_hat;
       for i = 1:N
            if abs(b(i)) <= 2*deniser_parameter && abs(b(i)) > deniser_parameter         
                x_hat(i) = sign(b(i))*(abs(abs(b(i))-deniser_parameter)+(abs(b(i))-deniser_parameter))/2;
            else
                if abs(b(i)) > 2*deniser_parameter && abs(b(i)) <= a*deniser_parameter
                x_hat(i) = ((a-1)*b(i)-sign(b(i))*a*deniser_parameter)/(a-2);
                else
                    if  abs(b(i)) > a*deniser_parameter
                        x_hat(i) = b(i);
                    else
                        x_hat(i) = 0;
                    end
                end
            end
       end    

      z = y - A*x_hat +  1/delta * z*1/N *sum(etaprimescad(b, deniser_parameter));
  end

     tau2(t) = 1/M*(norm(z,2)^2);
 
end


x_pre = x_hat;



