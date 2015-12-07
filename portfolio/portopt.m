load data_portfolio;
 
r = 100*r_cl(:,ind+1);
 
N = size(r,2);
T = size(r,1);
 
w_var1_lasso      = NaN*ones(455,N);
w_var1_postlasso  = NaN*ones(455,N);
w_var1_adalasso   = NaN*ones(455,N);
w_var1_lasso_logm = NaN*ones(455,N);
 
w_var20_lasso      = NaN*ones(455,N);
w_var20_postlasso  = NaN*ones(455,N);
w_var20_adalasso   = NaN*ones(455,N);
w_var20_lasso_logm = NaN*ones(455,N);
 
w_nochange_cens = NaN*ones(455,N);
w_dcc           = NaN*ones(455,N);
w_ewma          = NaN*ones(455,N);
 
t0   = 1023;
Ra   = 0.1;
mu_p = ((1+Ra).^(1/252)-1)*100;

for t = 1:455
    disp(t)
    
    
    %% VAR(1) LASSO
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = ivech(S_VAR1_LASSO(t,:)');
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_var1_lasso(t,:) = computewr(N,S,mu,mu_p); 
    % *******************************************************************
    
    
    %% VAR(1) POST LASSO
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = ivech(S_VAR1_POSTLASSO(t,:)');
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_var1_postlasso(t,:) = computewr(N,S,mu,mu_p); 
    % *******************************************************************
    
    %% VAR(1) ADALASSO
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = ivech(S_VAR1_ADALASSO(t,:)');
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_var1_adalasso(t,:) = computewr(N,S,mu,mu_p);
    % *******************************************************************
    
    %% VAR(1) LASSO LOGMAT
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = expm(ivech(S_VAR1_LASSO_LOGMAT(t,:)'));
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_var1_lasso_logm(t,:) = computewr(N,S,mu,mu_p); 
    % *******************************************************************
    
    %% VAR(20) LASSO
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = ivech(S_VAR20_LASSO(t,:)');
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_var20_lasso(t,:) = computewr(N,S,mu,mu_p); 
    % *******************************************************************
    
    
    %% VAR(20) POST LASSO
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = ivech(S_VAR20_POSTLASSO(t,:)');
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_var20_postlasso(t,:) = computewr(N,S,mu,mu_p); 
    % *******************************************************************
    
    %% VAR(20) ADALASSO
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = ivech(S_VAR20_ADALASSO(t,:)');
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_var20_adalasso(t,:) = computewr(N,S,mu,mu_p);
    % *******************************************************************
    
    %% VAR(20) LASSO LOGMAT
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = expm(ivech(S_VAR20_LASSO_LOGMAT(t,:)'));
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_var20_lasso_logm(t,:) = computewr(N,S,mu,mu_p); 
    % *******************************************************************
    
    %% NO-CHANGE CENSURED
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = ivech(S_NOCHANGE_CENS(t,:)');
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_nochange_cens(t,:) = computewr(N,S,mu,mu_p);
    % *******************************************************************
    
    %% DCC
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = ivech(S_DCC(t,:)');
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_dcc(t,:) = computewr(N,S,mu,mu_p);
    % *******************************************************************
    
    %% EWMA
    mu = mean(r(t0-2+t:-1:t0-2+t-100,:),1)';
    S  = ivech(S_EWMA(t,:)');
    
    [V,D]    = eig(S);
    dD       = diag(D);
    minD     = min(dD(dD>0));
    dD(dD<0) = minD;
    D        = diag(dD);
    S        = V*(D/V);
    
    w_ewma(t,:) = computewr(N,S,mu,mu_p);
    % *******************************************************************
end

w_lasso(:,:,1)  = w_var1_lasso;
w_lasso(:,:,2)  = w_var1_postlasso;
w_lasso(:,:,3)  = w_var1_adalasso;
w_lasso(:,:,4)  = w_var1_lasso_logm;
w_lasso(:,:,5)  = w_var20_lasso;
w_lasso(:,:,6)  = w_var20_postlasso;
w_lasso(:,:,7)  = w_var20_adalasso;
w_lasso(:,:,8)  = w_var20_lasso_logm;

r = r/100;
computeret
makestat
maketables

save results_day