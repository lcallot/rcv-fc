function [s,g] = mean_variance(w,S)

s = w'*S*w;
g = 2*w'*S;