function delta =  ecovalue(r1,r2,gamma)

delta = 0;

OptimizerOptions = optimset('Display','iter','GradObj','off','GradConstr','off','MaxFunEvals',1e10,...
                            'LargeScale','off','MaxIter',1000,'TolFun',1e-10,...
                            'DerivativeCheck','off','FunValCheck', 'on',...
                            'TolX',1e-10,'Diagnostics','off','FinDiffType','central',...
                            'UseParallel','always');
                        
delta = fsolve(@(delta)diff_util(delta,r1,r2,gamma),delta,OptimizerOptions);

delta = (1+delta).^(252)-1;
