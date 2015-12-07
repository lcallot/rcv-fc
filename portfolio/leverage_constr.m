function [c,ceq] = leverage_constr(w)

c   = sum(abs(w(w<0)))-0.3;
ceq = [];