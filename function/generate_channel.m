function [h] = generate_channel(para, path_loss)
%Generate the BS-user channels 
%  [h] = generate_channel(para, angle, path_loss)
%Inputs:
%   para: structure of the initial parameters
%   user_angle: struture of the angles
%   path_loss: structure of the path loss
%Outputs:
%   d: BS-user channels
%Date: 14/07/2021
%Edit: 27/09/2021
%Author: Zhaolin Wang

h = 1/sqrt(2) .* ( randn(para.N,para.K) + 1i*randn(para.N,para.K) );
h_normlized = h ./ vecnorm(h);

% correlation matrix
R = eye(para.K)/2;
phi = rand(1) * 2*pi;
for i = 1:para.K
    for j = i+1:para.K       
        R(i,j) = (para.t * exp(1i * phi))^(j-i);       
    end
end
R = R + R';

h_normlized = h_normlized * R^(1/2);
h = h_normlized .* vecnorm(h);

h = path_loss .* h;

%% re-order channel (from min to max)
h_norm = sum(abs(h).^2);
[~,I] = sort(h_norm);
h = h(:,I);

end

