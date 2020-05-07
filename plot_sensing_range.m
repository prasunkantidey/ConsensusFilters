function plot_sensing_range( x, y, r )
%PLOT_CIRCLE Summary of this function goes here
%   Detailed explanation goes here

    for i = 1:length(x)
        ang=0:0.01:2*pi; 
        xp=r(i)*cos(ang);
        yp=r(i)*sin(ang);
        plot(x(i)+xp,y(i)+yp);
    end
end
