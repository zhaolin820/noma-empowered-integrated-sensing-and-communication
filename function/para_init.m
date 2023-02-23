function [values] = para_init()
%Construct a struct of the initial values for all the parameters 
%  [values] = para_init()
%Inputs:
%   None
%Outputs:
%   values: a struct
%Date: 28/05/2021
%Author: Zhaolin Wang

values.noise = -120; % noise powe in dBm

values.N = 8; % overall antennas

values.Pt = 10^(20/10); % overall transmit power
values.n = 1; % equivalent noise power
values.K = 3; % user number
values.theta = [-40, 40]*pi/180; % desired directions
values.weight = ones(values.K,1); % weight of weighted sum rate 
values.pathloss_direct =  @(d) 32.6 + 36.7*log10(d); % path loss with d in m

values.BS_loc = [0,0]; % location of BS
values.user_center = [0, 0];
values.user_range = [50, 200]; % range of user location from RIS in m

values.t = 0; % channel correlation factor
values.R_min = 1; % minimum rate
values.correlation_min = 10; % minimum cross-correlation

end




