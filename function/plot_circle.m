function plot_circle(o,r)
%Plot a circle
%   plotCircle(o,r)
%Inputs:
%   o: Center of circle
%   r: Radius of circle
%Outputs:
%   None
%Date: 27/12/2020
%Author: Zhaolin Wang

angle = 0:0.001:2*pi;
x = r*cos(angle) + o(1);
y = r*sin(angle) + o(2);
c = 0.7*[1 1 1];
plot(x,y,'Color',c);
end

