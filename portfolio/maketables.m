
avg_W    = [squeeze(mean(mean(abs(w_lasso))))' squeeze(mean(mean(abs(w_nochange_cens))))' squeeze(mean(mean(abs(w_dcc))))' squeeze(mean(mean(abs(w_ewma))))'];
max_W    = [squeeze(max(max((w_lasso))))' squeeze(max(max((w_nochange_cens))))' squeeze(max(max((w_dcc))))' squeeze(max(max((w_ewma))))'];
min_W    = [squeeze(min(min((w_lasso))))' squeeze(min(min((w_nochange_cens))))' squeeze(min(min((w_dcc))))' squeeze(min(min((w_ewma))))'];
avg_lev  = [mean(squeeze(sum(w_lasso.*(w_lasso<0),2))) mean(squeeze(sum(w_nochange_cens.*(w_nochange_cens<0),2))) mean(squeeze(sum(w_dcc.*(w_dcc<0),2)))  mean(squeeze(sum(w_ewma.*(w_ewma<0),2)))];
prop_lev = [mean(squeeze(mean((w_lasso<0),2))) mean(squeeze(mean((w_nochange_cens<0),2))) mean(squeeze(mean((w_dcc<0),2))) mean(squeeze(mean((w_ewma<0),2)))];
avg_turn = [mean(squeeze(mean(abs(w_lasso-w_lasso_hold),2)))  mean(squeeze(mean(abs(w_nochange_cens-w_nochange_hold),2))) mean(squeeze(mean(abs(w_dcc-w_dcc_hold),2)))  mean(squeeze(mean(abs(w_ewma-w_ewma_hold),2)))];

aux = NaN*ones(2,3);
TABLE = [avg_W;max_W;min_W;avg_lev;prop_lev;avg_turn;mu_p_R;cumr_p_R;std_p_R;sr_p_R;dr_p_R;100*[squeeze(Delta(:,:,1)') aux];100*[squeeze(Delta(:,:,2)') aux];100*[squeeze(Delta(:,:,3)') aux]]