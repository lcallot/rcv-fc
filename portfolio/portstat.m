function [mu_p,cumr_p,std_p,var_p,sr_p,dr_p] =  portstat(r,w,vS)

T = size(r,1);

dr_p = NaN*ones(T,1);

mu_p   = mean(r);
cumr_p = prod(1+r);
std_p  = std(r);
var_p  = var(r);
sr_p   = mu_p/std_p;

for t = 1:T
    S        = ivech(vS(t,:)');
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    aux = diag(S);
    s   = sqrt(w(t,:)*S*w(t,:)');
    dr_p(t,1) = (w(t,:)*sqrt(aux))/s;
end

dr_p = mean(dr_p);