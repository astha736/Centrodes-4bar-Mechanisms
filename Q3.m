% Author: Astha Gupta 
% Student Id: S4899512

clc;
clear all;
close all;

% Parameters given 
n=6;  
la = 4;
lb = 1/(0.1+(0.012*n));

% fixed points
Oa = [0,0]';
A = [0,lb]';

figure(1)

for i=1:10
% calculating the Contact point as the mechanism moves 
for alpha_p= 0:0.1:2*pi
    
alpha = pi/2 - alpha_p;
% Calculating theta and beta angle given alpha
Ob_A = sqrt((la^2)+(lb^2) -(2*la*lb*(cos(alpha))));
ck = (la*sin(alpha))/Ob_A;
kangle = asin(ck); % in triangle Ob-B-P the Ap = ObP => isosceles trainagle 
% angle  Ob-B-P + theta = 180; angle on the same line segment OaA
% angle  Ob-B-P + 2*k = 180; property of trainagle 
theta = 2*kangle
beta = (pi-theta-alpha);


% calculate the Point P (IC and point of contact of the two ellipse)
m1 = (sin(beta)*la)/sin(theta); %sine rule and Oa-P-Ob traingle 
% P = m1*[cos(alpha),sin(alpha)]';
P = [0,m1];

% A = [lb*cos(alpha) lb*sin(alpha)];
Ob = [la*cos(alpha_p) la*sin(alpha_p)];
B = [lb*cos(pi/2+theta)+Ob(1) lb*sin(pi/2+theta)+Ob(2)];

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

%  contact point at intersection of line through AB, OaOb 
x1= [A(1) B(1)];
y1= [A(2) B(2)];

x2= [Oa(1) Ob(1)];
y2= [Oa(2) Ob(2)];

ply1 = polyfit(x1,y1,1);
ply2 = polyfit(x2,y2,1);

x_intersect = fzero(@(x) polyval(ply1-ply2,x),3);
y_intersect = polyval(ply1,x_intersect);
plot(x_intersect,y_intersect,'r*');
hold on

X_mov = zeros(63,1);
Y_mov = zeros(63,1);
i = 0;

    % calculating the plot for moving centrode for a given configuration!
    for alpha_k= 0:0.1:2*pi

    alpha_n = alpha + alpha_k;
    % Calculating theta and beta angle given alpha
    Ob_A_n = sqrt((la^2)+(lb^2) -(2*la*lb*(cos(alpha_n))));
    ck_n = (la*sin(alpha_n))/Ob_A_n;
    kangle_n = asin(ck_n); % in triangle Ob-B-P the Ap = ObP => isosceles trainagle 
    % angle  Ob-B-P + theta = 180; angle on the same line segment OaA
    % angle  Ob-B-P + 2*k = 180; property of trainagle 
    theta_n = 2*kangle_n;
    beta_n = (pi-theta_n-alpha_n);

    beta_phi = beta - alpha_p;
    A_n = [la*cos(alpha_n - beta_phi)+B(1) la*sin(alpha_n - beta_phi)+B(2) ];
    Oa_n = [la*cos(pi+beta_n - beta_phi)+Ob(1) la*sin(pi+beta_n - beta_phi)+Ob(2)];

    x1= [A_n(1) B(1)];
    y1= [A_n(2) B(2)];

    x2= [Oa_n(1) Ob(1)];
    y2= [Oa_n(2) Ob(2)];

    ply1 = polyfit(x1,y1,1);
    ply2 = polyfit(x2,y2,1);
    i = i+1;
    X_mov(i) = fzero(@(x) polyval(ply1-ply2,x),3);
    Y_mov(i) = polyval(ply1,X_mov(i));

    end

p5 = plot(X_mov,Y_mov,'-o');

getframe;
% 
delete(p1);
delete(p2);
delete(p3);
delete(p4);
delete(p5);
% delete(p5);
% delete(p6);


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
% p5 = plot(xe_mov, ye_mov,'b');

