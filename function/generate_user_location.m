function [user_loc, user_angle, d_BU] = generate_user_location(para)
%Generate the user locations randomly 
%  [user_loc, user_angle, d_RU, d_BU] = generate_user_location(para)
%Inputs:
%   para: structure of the initial parameters
%Outputs:
%   user_loc: locations of users in meters
%   user_agnle: angle of users from RIS
%   d_BU: distance between BS and RIS
%Date: 30/05/2021
%Edit: 27/09/2021
%Author: Zhaolin Wang


user_gap = (para.user_range(2) - para.user_range(1)) / (para.K-1);
d_CU = (0:para.K-1)' * user_gap + para.user_range(1);
user_angle_center = rand(para.K,1) * (2 *pi); % angle of directions from center

user_loc = d_CU.*exp(1i*user_angle_center);
user_loc = [real(user_loc), imag(user_loc)] + para.user_center;

relative_user_loc = user_loc - para.BS_loc;
d_BU = sqrt(relative_user_loc(:,1).^2 + relative_user_loc(:,2).^2);
user_angle = atan(relative_user_loc(:,2) ./ relative_user_loc(:,1));



end

