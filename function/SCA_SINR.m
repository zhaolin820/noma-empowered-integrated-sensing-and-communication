function [sca_all] = SCA_SINR(para, H, W, W_n, r)
%The SINR constraints in SCA. 
%  [sca] = SCA_SINR(para, H, W, W_n, r, r_n)
%Inputs:
%   para: system paramaters
%   H: channel gain
%   W: beamformer (cvx)
%   W_n beamformer in the last step
%   r: auxiliary variable (cvx)
%   r_n: auxiliary variable in the last step 
%Outputs:
%   sca: convex cvx variables for SINR constraints
%Date: 06/10/2021
%Author: Zhaolin Wang

sca_all = [];
for k = 1:para.K-1
    W_k = W(:,:,k);
    r_k = r(k);
    
    % decoding
    for j = k:para.K
        H_j = H(:,:,j);
        
        % interference
        I = 0; I_n = 0; I_diff = 0;
        for i = k+1:para.K
            W_i = W(:,:,i);
            W_n_i = W_n(:,:,i);
            I = I + real(trace(H_j*W_i));
            I_n = I_n +  real(trace(H_j*W_n_i));
            I_diff = I_diff + real(trace(H_j * (W_i - W_n_i)));
        end
        F_j_k = - log(para.n + I_n)/log(2) - I_diff / (para.n + I_n) / log(2);
        sca = log(para.n + real(trace(H_j * W_k)) + I)/log(2) + F_j_k - r_k;
        sca_all = [sca_all, sca];
    end   
end

H_K = H(:,:,end); W_K = W(:,:,end); r_K = r(end);
sca_all = [sca_all, log( 1 + real(trace(H_K*W_K))/para.n )/log(2)-r_K];
end
