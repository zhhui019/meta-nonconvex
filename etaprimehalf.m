function y = etaprimehalf(b, lambda_t)
% etaprimehalf is the derivative of the half threshold function. 

N = length(b);
y = zeros(N,1);


 

 for i = 1:N
      if abs(b(i)) >= 3/4*(lambda_t^(2/3))
          phi = acos(lambda_t/8*((abs(b(i))/3)^(-1.5)));
          y(i) =4*b(i)/3*(cos(((pi)/3)-phi/3))^2+(((sqrt(3))/4)*lambda_t*(sin(2/3*(pi - phi))))/sqrt(b(i)^3-(27/64)*lambda_t^2);
      else
          y(i) = 0;
      end
 end
             