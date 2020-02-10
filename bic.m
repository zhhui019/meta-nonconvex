% select lambda by BIC


function value = bic(X,Y,obs,coef)

start_idx = cumsum(obs) + 1 - obs;
end_idx  = cumsum(obs) ;            
  M      = size(coef,2);
  value  = 0;
  
  for m = 1:M
      prob = X(start_idx(m):end_idx(m),:)*coef(:,m);
      prob = exp(prob)./(1+exp(prob));
      for i=1:length(prob)
          if prob(i)<1e-10
              prob(i)=1e-10;
          elseif prob(i)>1-(1e-10)
              prob(i)=1-(1e-10);
          end
      end
      p=sum(coef(:,m)~=0);
      s_tmp = Y(start_idx(m):end_idx(m)).*log(prob)+(1-Y(start_idx(m):end_idx(m)).*log(1-prob));
      value = value+sum(s_tmp)/obs(m)*(-2)+p*(log(obs(m))+2*log(size(X,2)))/obs(m);
      
  end