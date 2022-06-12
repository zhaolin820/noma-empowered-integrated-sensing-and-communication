function plot_location(para, user_loc)
%Plot the locations of BS, RIS and users
%  plot_location(para, user_loc)
%Inputs:
%   para: structure of the initial parameters
%   user_loc: locations of users
%Outputs:
%   None
%Date: 30/05/2021
%Author: Zhaolin Wang

figure;
plot(para.BS_loc(1), para.BS_loc(2), '^r','MarkerSize',10,'LineWidth',2);
hold on;
plot(user_loc(:,1), user_loc(:,2),'ok','MarkerSize',6,'LineWidth',2);
plot_circle(para.user_center,para.user_range(1));
plot_circle(para.user_center,para.user_range(2))
xlim([-200,200]);
ylim([-200,200]);
axis equal;
grid on;


end

