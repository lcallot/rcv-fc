function w = computewr(N,S,r,mu_p)

lb = -0.20*ones(N,1);
ub = 0.20*ones(N,1);

Aeq = [ones(1,N);r'];
beq = [1;mu_p];
A   = [];
b   = [];

w = ones(N,1)*(1/N);

OptimizerOptions = optimset('Display','notify','GradObj','on','GradConstr','off','MaxFunEvals',1e8,...
                            'LargeScale','off','MaxIter',2000,'TolFun',1e-5,...
                            'DerivativeCheck','off','FunValCheck', 'on',...
                            'TolX',1e-6,'Algorithm','interior-point','Diagnostics','off','FinDiffType','central',...
                            'UseParallel','always','SubproblemAlgorithm','cg');

w = fmincon(@(w)mean_variance(w,S),w,A,b,Aeq,beq,lb,ub,@leverage_constr,OptimizerOptions);
