% Bowling Ball Dynamics with Euler Angles
% DAE Solution

clear
close all
clc

% FUNCTION HANDLES FOR SYMBOLIC SOLUTIONS
fsymbolicSlip = slippingDAE();
fsymbolicRoll = rollingDAE();

%%
% % PARAMS REGULAR BALL
% p.g = 9.81; % m s^-2
% p.R = .1080; % m
% p.m = 6.35; % kg
% p.I1 = .001; % kg m^2
% p.I2 = .001; % kg m^2
% p.I3 = .03; % kg m^2
% 
% % INITS REGULAR BALL
% x0 = 0; y0 = 0; z0 = p.R; % m
% x_dot0 = 8; y_dot0 = -.45; z_dot0 = 0;% m/s
% phi0 = deg2rad(90); theta0 = deg2rad(90); psi0 = deg2rad(0); % rad
% phi_dot0 = 0; theta_dot0 = 0; psi_dot0 = -30; % rad/s
% position0 = [x0 y0 z0 x_dot0 y_dot0 z_dot0]';
% attitude0 = [phi0 theta0 psi0 phi_dot0 theta_dot0 psi_dot0]';
% X0 = [position0;attitude0];

%%
% PARAMS BIGGER BALL
p.g = 9.81; % m s^-2
p.R = .3080; % m
p.m = 6.35; % kg
p.I1 = .03;
p.I2 = .03;
p.I3 = .2;

% INITS BIGGER BALL
x0 = 0; y0 = 0; z0 = p.R; % m
x_dot0 = 8; y_dot0 = -.6; z_dot0 = 0;% m/s
phi0 = deg2rad(90); theta0 = deg2rad(90); psi0 = deg2rad(0); % rad
phi_dot0 = 0; theta_dot0 = 0; psi_dot0 = -16; % rad/s
position0 = [x0 y0 z0 x_dot0 y_dot0 z_dot0]';
attitude0 = [phi0 theta0 psi0 phi_dot0 theta_dot0 psi_dot0]';
X0 = [position0;attitude0];

%% TIME
step = 0.01; tend = 2.5; % s
tode = 0:step:tend;

% ODE45
f = @(tode,X) bowlingballDAE(tode,X,p,fsymbolicSlip,fsymbolicRoll);
options = odeset('absTol',1e-6,'relTol',1e-6);
[t,X] = ode45(f,tode,X0);

% CALCULATE REACTIONS
len = length(X);
phi = X(:,7); theta = X(:,8);
phi_dot = X(:,10); theta_dot = X(:,11); psi_dot = X(:,12);
x_dot = X(:,4); y_dot = X(:,5); z_dot = X(:,6);
Nrc = [0;0;-p.R];
vco = zeros(3,len);
vconorm = zeros(len,1);
for i = 1:len
    omega_x = theta_dot(i)*cos(phi(i))+psi_dot(i)*sin(theta(i))*sin(phi(i));
    omega_y = theta_dot(i)*sin(phi(i))-psi_dot(i)*sin(theta(i))*cos(phi(i));
    omega_z = phi_dot(i)+psi_dot(i)*cos(theta(i));
    Nomega = [omega_x;omega_y;omega_z];
    v = [x_dot(i);y_dot(i);z_dot(i)];
    vcg = cross(Nomega,Nrc);
    vgo = v;
    vco(:,i) = vcg+vgo;
    vconorm(i) = norm(vco(:,i));
end

%{
plot(t,vconorm,'k-',LineWidth=1);
xlabel("t (s)", Interpreter="latex");
ylabel("$v_{C/O} (\frac{m}{s})$", Interpreter="latex");
title("$v_{C/O}$ vs time",Interpreter="latex")
set(gca,'fontsize',16, 'fontname','palatino linotype'); grid;
%}

% ANIMATE
mode = 1; % 1 = animation, 2 = animation + record & save it
hgt_animate(X,step,p,mode)

function Xdot = bowlingballDAE(t,X,p,fsymbolicSlip,fsymbolicRoll)
% X = [x y z x_dot y_dot z_dot
%      phi theta psi phi_dot theta_dot psi_dot]'
% Xdot = [x_dot y_dot z_dot x_ddot y_ddot z_ddot
%         phi_dot theta_dot psi_dot phi_ddot theta_ddot psi_ddot]'

% PARAMS
g = p.g; R = p.R; m = p.m;
I1 = p.I1; I2 = p.I2; I3 = p.I3;
Nrc = [0;0;-R];

% COORDS
x = X(1);
x_dot = X(4); y_dot = X(5); z_dot = X(6);
v = [x_dot;y_dot;z_dot];
phi = X(7); theta = X(8);
phi_dot = X(10); theta_dot = X(11); psi_dot = X(12);

% DETERMINE FRICTION
length = 18.288;
if x < .5*length
    mu = .04;
elseif x >= .5*length && x <= .7*length
    mu = 1.3/length*x-.61;
else
    mu = .3;
end

% EULER ANGLE RATES --> INERTIAL AXIS RATES
omega_x = theta_dot*cos(phi)+psi_dot*sin(theta)*sin(phi);
omega_y = theta_dot*sin(phi)-psi_dot*sin(theta)*cos(phi);
omega_z = phi_dot+psi_dot*cos(theta);
Nomega = [omega_x;omega_y;omega_z];
vcg = cross(Nomega,Nrc);
vgo = v;
vco = vcg+vgo; % velocity of contact point relative to origin
if norm(vco) > .001 % slipping
    Xdot = fsymbolicSlip(t,X,vco(1),vco(2),vco(3),g,m,R,mu,I1,I2,I3);
else % rolling
    Xdot = fsymbolicRoll(t,X,g,m,R,mu,I1,I2,I3);
end

end
