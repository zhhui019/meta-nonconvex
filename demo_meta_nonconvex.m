%the simulation demo for meta-analysis method based on nonconvex regularization

clear all
n = 50;
p = 1000;
M = 10;
obs = repelem(n,M);
maxit = 30;
X = [];
Y = [];
beta_true = [];

for m=1:M
        X_tmp = randn(n,p);  
        beta =  [1, -1, 2, -1, zeros(1,p - 4)]';
        pb = X_tmp*beta;
        prob = exp(pb)./(1+exp(pb));
        Y_tmp = zeros(n,1);
		Y_tmp = prob>0.5;
        X =[X;X_tmp];
        Y = [Y;Y_tmp];
        beta_true = [beta_true;beta];
 
 end
 
 %select lambda
lams = 0.01:0.01:0.08;
BIC=zeros(1,length(lams));
for j = 1:length(lams)
mycoef_bic = meta_nonconvex(X,Y, obs,lams(j), maxit,'scad',1e-3);
BIC(j) = bic(X,Y,obs,mycoef_bic);    
end
[min_lam, idx_min] = min(BIC);
bestfit_coef = meta_nonconvex(X,Y, obs,lams(idx_min), maxit,'scad',1e-3);