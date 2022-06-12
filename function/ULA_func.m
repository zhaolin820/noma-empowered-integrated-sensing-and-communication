function [a] = ULA_func(theta, N)
%The steering vector of ULA with half wavelength distancing
%  [a] = ULA_func(theta, N)
%Inputs:
%   theta: angle of direction
%   N: number of antennas
%Outputs:
%   a: steering vector
%Date: 28/02/2021
%Author: Zhaolin Wang

n = 0:(N-1);
a = exp(1i * pi * n * sin(theta)).';

end

