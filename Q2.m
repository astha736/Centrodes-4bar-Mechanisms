% Author: Astha Gupta 
% Student Id: S4899512

clc;
clear all;
close all;

% Parameters given 
n=6;  
la = 4;
lb = 1/(0.1+(0.012*n));

Oa = [0,0]';
Ob = [la,0]';

figure(1)

% Drawing our fixed centrodes
Cx = (Oa(1)+Ob(1))/2;
Cy = (Oa(2)+Ob(2))/2;
r1 = lb/2; % semi-axis length 
r2 = sqrt(lb^2 - la^2)/2; % semi-axis length 
angle = 0;
[xe_fixed,ye_fixed] = ellipseParam(Cx,Cy,r1,r2,angle);

pF = plot(xe_fixed, ye_fixed,'g');
hold on

for i=1:10
% calculating the Contact point as the mechanism moves 
for alpha= 0:0.1:2*pi

% Calculating theta and beta angle given alpha
Ob_B = sqrt((la^2)+(lb^2) -(2*la*lb*(cos(alpha))));
ck = (la*sin(alpha))/Ob_B;
kangle = asin(ck); % in triangle Ob-B-P the Ap = ObP => isosceles trainagle 
% angle  Ob-B-P + theta = 180; angle on the same line segment OaA
% angle  Ob-B-P + 2*k = 180; property of trainagle 
theta = 2*kangle
beta = (pi-theta-alpha);


% calculate the Point P (IC and point of contact of the two ellipse)
m1 = (sin(beta)*la)/sin(theta); %sine rule and Oa-P-Ob traingle 
P = m1*[cos(alpha),sin(alpha)]';

A = [lb*cos(alpha) lb*sin(alpha)];
B = [la-lb*cos(beta) lb*sin(beta)]; % cos(180-beta) = -cos(beta); sin(180-beta) = -sin(beta) 

% Calculting our moving centrode 
Cx = (A(1)+B(1))/2;
Cy = (A(2)+B(2))/2;
r1 = lb/2; % semi-axis length 
r2 = sqrt(lb^2 - la^2)/2; % semi-axis length 
slope = (B(2)- A(2))/ (B(1)-A(1));
angle = atan(slope);

[xe_mov,ye_mov] = ellipseParam(Cx,Cy,r1,r2,angle);


% x-value and v-values 
axis([-10 10 -10 10]);
% plotting the mechanism
p1 = plot([Oa(1) A(1)], [Oa(2) A(2)]);
hold on
p2 = plot([Oa(1) Ob(1)], [Oa(2) Ob(2)]);
hold on
p3 = plot([Ob(1) B(1)], [Ob(2) B(2)]);
hold on
p4 = plot([A(1) B(1)], [A(2) B(2)]);
hold on

% contact point P
plot(P(1),P(2),'r*');
hold on

% moving centrode
p5 = plot(xe_mov, ye_mov,'b');
hold on

% Center of rod AB
xc = (A(1)+B(1))/2;
yc = (A(2)+B(2))/2;
p6 = plot(xc,yc,'b*');

getframe;

delete(p1);
delete(p2);
delete(p3);
delete(p4);
delete(p5);
delete(p6);


end

end

p1 = plot([Oa(1) A(1)], [Oa(2) A(2)]);
hold on
p2 = plot([Oa(1) Ob(1)], [Oa(2) Ob(2)]);
hold on
p3 = plot([Ob(1) B(1)], [Ob(2) B(2)]);
hold on
p4 = plot([A(1) B(1)], [A(2) B(2)]);
hold on 
p5 = plot(xe_mov, ye_mov,'b');

