function [P] = probing_power(para, Rw)
%Calculate the probing power (effective sensing power)
%  [P] = probing_power(para, Rw)
%Inputs:
%   para: system paramaters
%   Rw: transmit covariance matrix
%Outputs:
%   P: probing power
%Date: 03/10/2021
%Author: Zhaolin Wang

P = [];
for m = 1:length(para.theta)
    a_m = ULA_func(para.theta(m), para.N);
    P = [P, real(a_m' * Rw * a_m)]; 
end
end

