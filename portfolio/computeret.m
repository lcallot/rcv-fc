w_lasso_hold    = NaN*ones(455,N,8);
w_nochange_hold = NaN*ones(455,N);
w_dcc_hold      = NaN*ones(455,N);
w_ewma_hold     = NaN*ones(455,N);

r_p_lasso    = NaN*ones(455,8);
r_p_nochange = NaN*ones(455,1);
r_p_dcc      = NaN*ones(455,1);
r_p_ewma     = NaN*ones(455,1);

omega0 = (1/30)*ones(1,30);

c = 0.001*ones(455,1);

for t = 1:455
    for i = 1:8
        if t == 1
            w_lasso_hold(t,:,i)    = omega0.*((1+r(t0+t-2,:))./(1+mean(r(t0+t-2,:))));
        else
            w_lasso_hold(t,:,i)    = w_lasso(t-1,:,i).*((1+r(t0+t-2,:))./(1+r_p_lasso(t-1,i)));
        end
    
        r_p_lasso(t,i)    = w_lasso(t,:,i)*(r(t0+t-1,:) - squeeze(c(t)*abs(w_lasso(t,:,i)-w_lasso_hold(t,:,i))))';
    end
    
    if t == 1
        w_nochange_hold(t,:) = omega0.*((1+r(t0+t-2,:))./(1+mean(r(t0+t-2,:))));
        w_dcc_hold(t,:)      = omega0.*((1+r(t0+t-2,:))./(1+mean(r(t0+t-2,:))));
        w_ewma_hold(t,:)     = omega0.*((1+r(t0+t-2,:))./(1+mean(r(t0+t-2,:))));
    else
        w_nochange_hold(t,:) = w_nochange_cens(t-1,:).*((1+r(t0+t-2,:))./(1+r_p_nochange(t-1)));
        w_dcc_hold(t,:)      = w_dcc(t-1,:).*((1+r(t0+t-2,:))./(1+r_p_dcc(t-1)));
        w_ewma_hold(t,:)     = w_ewma(t-1,:).*((1+r(t0+t-2,:))./(1+r_p_ewma(t-1)));
    end
    r_p_nochange(t,1) = w_nochange_cens(t,:)*(r(t0+t-1,:) - squeeze(c(t)*abs(w_nochange_cens(t,:)-w_nochange_hold(t,:))))';
    r_p_dcc(t,1)      = w_dcc(t,:)*(r(t0+t-1,:) - squeeze(c(t)*abs(w_dcc(t,:)-w_dcc_hold(t,:))))';   
    r_p_ewma(t,1)     = w_ewma(t,:)*(r(t0+t-1,:) - squeeze(c(t)*abs(w_ewma(t,:)-w_ewma_hold(t,:))))';   
end    
