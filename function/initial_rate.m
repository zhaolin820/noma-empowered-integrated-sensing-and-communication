function [ineqn] = initial_rate(para,H, W)
%For rate initialization in SCA algorithm
%  [ineqn] = initial_rate(para,H, W)
%Inputs:
%   para: system paramaters
%   H: channel gains
%   W: beamformer
%Outputs:
%   ineqn: cvx inequalities
%Date: 03/10/2021
%Author: Zhaolin Wang

ineqn = [];
for k = 1:para.K-1
    W_k = W(:,:,k);
    for j = k:para.K
        H_j = H(:,:,j);
        
        % interference
        I = 0;
        for i = k+1:para.K
            W_i = W(:,:,i);
            I = I + real(trace(H_j*W_i));
        end
        ineqn = [ineqn, -real(trace(H_j*W_k)) + (2^para.R_min -1)*(I + para.n)];
    end   
end

H_K = H(:,:,end); W_K = W(:,:,end);
ineqn = [ineqn, -real(trace(H_K*W_K)) + (2^para.R_min - 1) * para.n];
end

