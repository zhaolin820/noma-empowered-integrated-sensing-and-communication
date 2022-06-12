function [corr] = cross_correlation(para, Rw)
%Calculate the cross correlation in different target directions
%  [P] = probing_power(para, Rw)
%Inputs:
%   para: system paramaters
%   Rw: transmit covariance matrix
%Outputs:
%   P: probing power
%Date: 03/10/2021
%Author: Zhaolin Wang


M = length(para.theta);
corr = 0;
for k = 1:M
    for p = k+1:M
        a_k = ULA_func(para.theta(k), para.N);
        a_p = ULA_func(para.theta(p), para.N);
        corr = corr + pow_abs(a_k' * Rw * a_p, 2);
    end
end
end

