function [R] = rate_NOMA(para, W, H)
%Calculate the achievable rate of NOMA-DFRC system
%  [R] = rate_NOMA(para, W, H)
%Inputs:
%   para: system paramaters
%   W: beamformer
%   H: channel gains
%Outputs:
%   R: achievable rate
%Date: 03/10/2021
%Author: Zhaolin Wang

R = inf * ones(para.K, para.K);

for k = 1:para.K-1
   W_k = W(:,:,k);
   
   % decoding user-k's stream at user-j
   for j = k:para.K
       H_j = H(:,:,j);
       
       % interference
       I = 0;
       for i = k+1:para.K
           W_i = W(:,:,i);
           I = I + real(trace(H_j * W_i));
       end     
       R(k,j) = log2( 1 + real(trace(H_j * W_k)) / (I + para.n) );
   end
end
W_K = W(:,:,end); H_K = H(:,:,end);
R(end,end) = log2( 1 + real(trace(H_K * W_K)) );
end

