% Descriptive Statistics

[mu_p_R(1),cumr_p_R(1),std_p_R(1),var_p_R(1),sr_p_R(1),dr_p_R(1)] = portstat(r_p_lasso(:,1),w_lasso(:,:,1),S_TRUE);
[mu_p_R(2),cumr_p_R(2),std_p_R(2),var_p_R(2),sr_p_R(2),dr_p_R(2)] = portstat(r_p_lasso(:,2),w_lasso(:,:,2),S_TRUE);
[mu_p_R(3),cumr_p_R(3),std_p_R(3),var_p_R(3),sr_p_R(3),dr_p_R(3)] = portstat(r_p_lasso(:,3),w_lasso(:,:,3),S_TRUE);
[mu_p_R(4),cumr_p_R(4),std_p_R(4),var_p_R(4),sr_p_R(4),dr_p_R(4)] = portstat(r_p_lasso(:,4),w_lasso(:,:,4),S_TRUE);
[mu_p_R(5),cumr_p_R(5),std_p_R(5),var_p_R(5),sr_p_R(5),dr_p_R(5)] = portstat(r_p_lasso(:,5),w_lasso(:,:,5),S_TRUE);
[mu_p_R(6),cumr_p_R(6),std_p_R(6),var_p_R(6),sr_p_R(6),dr_p_R(6)] = portstat(r_p_lasso(:,6),w_lasso(:,:,6),S_TRUE);
[mu_p_R(7),cumr_p_R(7),std_p_R(7),var_p_R(7),sr_p_R(7),dr_p_R(7)] = portstat(r_p_lasso(:,7),w_lasso(:,:,7),S_TRUE);
[mu_p_R(8),cumr_p_R(8),std_p_R(8),var_p_R(8),sr_p_R(8),dr_p_R(8)] = portstat(r_p_lasso(:,8),w_lasso(:,:,8),S_TRUE);

[mu_p_R(9),cumr_p_R(9),std_p_R(9),var_p_R(9),sr_p_R(9),dr_p_R(9)]       = portstat(r_p_nochange,w_nochange_cens,S_TRUE);
[mu_p_R(10),cumr_p_R(10),std_p_R(10),var_p_R(10),sr_p_R(10),dr_p_R(10)] = portstat(r_p_dcc,w_dcc,S_TRUE);
[mu_p_R(11),cumr_p_R(11),std_p_R(11),var_p_R(11),sr_p_R(11),dr_p_R(11)] = portstat(r_p_ewma,w_ewma,S_TRUE);

gamma = [1;5;10];

Delta = NaN*ones(8,3,3);
for j = 1:length(gamma)
    r1 = r_p_nochange;
    for k = 1:8      
        r2           = r_p_lasso(:,k);
        Delta(k,1,j) = ecovalue(r1,r2,gamma(j));
    end
    
    r1 = r_p_dcc;
    for k = 1:8      
        r2           = r_p_lasso(:,k);
        Delta(k,2,j) = ecovalue(r1,r2,gamma(j));
    end
    
    r1 = r_p_ewma;
    for k = 1:8      
        r2           = r_p_lasso(:,k);
        Delta(k,3,j) = ecovalue(r1,r2,gamma(j));
    end
end

