function f = diff_util(delta,r1,r2,gamma)

f = mean((1+r1)-(gamma/(2*(1+gamma)))*(1+r1).^2) - mean((1+r2-delta)-(gamma/(2*(1+gamma)))*(1+r2-delta).^2);
