function [xr,yr] = ellipseParam(xc,yc,r1,r2,angle)
% Author: Astha Gupta 
% Student Id: S4899512

t = linspace(0, 2*pi, 200);
xt = r1 * cos(t) + xc;
yt = r2 * sin(t) + yc;
% % aply rotation by angle theta

R = [cos(angle) -sin(angle) ; sin(angle) cos(angle)] ;
P = zeros(2,length(xt)) ;
for i =1:length(xt)
    P(:,i) = R*[xt(i)-xc ;yt(i)-yc] + [xc;yc];
end
xr = P(1,:) ; yr = P(2,:) ;
end