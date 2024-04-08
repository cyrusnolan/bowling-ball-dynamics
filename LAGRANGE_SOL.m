% Bowling Ball Dynamics with Euler Angles
% Euler-Lagrange Solution

clear
close all
clc

% PARAMS - SLIPPING DEMO
p.g = 9.81; % m s^-2
p.R = .1080; % m
p.m = 6.35; % kg
p.I1 = .001; % kg m^2
p.I2 = .001; % kg m^2
p.I3 = .03; % kg m^2

% INITS - SLIPPING DEMO
x0 = 0; y0 = 0; z0 = p.R; % m
x_dot0 = 8; y_dot0 = -.5; z_dot0 = 0;% m/s
phi0 = deg2rad(90); theta0 = deg2rad(90); psi0 = deg2rad(0); % rad
phi_dot0 = 0; theta_dot0 = 0; psi_dot0 = -30; % rad/s
position0 = [x0 y0 z0 x_dot0 y_dot0 z_dot0]';
attitude0 = [phi0 theta0 psi0 phi_dot0 theta_dot0 psi_dot0]';
X0 = [position0;attitude0];
fsymbolic = slippingEL();

%% PARAmS BIGGER BALL
% p.g = 9.81; % m s^-2
% p.R = .3080; % m
% p.m = 6.35; % kg
% p.I1 = .03;
% p.I2 = .03;
% p.I3 = .2;
% 
% % INITS BIGGER BALL
% x0 = 0; y0 = 0; z0 = p.R; % m
% x_dot0 = 8; y_dot0 = -.6; z_dot0 = 0;% m/s
% phi0 = deg2rad(90); theta0 = deg2rad(90); psi0 = deg2rad(0); % rad
% phi_dot0 = 0; theta_dot0 = 0; psi_dot0 = -16; % rad/s
% position0 = [x0 y0 z0 x_dot0 y_dot0 z_dot0]';
% attitude0 = [phi0 theta0 psi0 phi_dot0 theta_dot0 psi_dot0]';
% X0 = [position0;attitude0];
% fsymbolic = slippingEL();

%% PARAMS - ROLLING DEMO
% p.g = 9.81; % m s^-2
% p.R = .3080; % m
% p.m = 6.35; % kg
% p.I1 = .05; % kg m^2
% p.I2 = .05; % kg m^2
% p.I3 = .05; % kg m^2
% 
% % INITS - ROLLING DEMO
% x0 = 0; y0 = .4; z0 = p.R; % m
% x_dot0 = 3; y_dot0 = -.2; z_dot0 = 0;% m/s
% phi0 = deg2rad(50); theta0 = deg2rad(90); psi0 = deg2rad(0); % rad
% wQe = [0 cos(phi0) sin(theta0)*sin(phi0);
%        0 sin(phi0) -sin(theta0)*cos(phi0);
%        1 0         cos(theta0)];
% omega_d = [-y_dot0/p.R;x_dot0/p.R;0];
% e_dot = wQe\omega_d;
% phi_dot0 = e_dot(1); theta_dot0 = e_dot(2); psi_dot0 = e_dot(3); % rad/s
% position0 = [x0 y0 z0 x_dot0 y_dot0 z_dot0]';
% attitude0 = [phi0 theta0 psi0 phi_dot0 theta_dot0 psi_dot0]';
% X0 = [position0;attitude0];
% fsymbolic = rollingEL();

%% TIME
step = 0.01; tend = 2.5; % s
tode = 0:step:tend;

% ODE45
f = @(tode,X) bowlingballEL(tode,X,p,fsymbolic);
options = odeset('absTol',1e-6,'relTol',1e-6);
[t,X] = ode45(f,tode,X0);

% CALCULATE REACTIONS
len = length(X);
x_dot = X(:,4); y_dot = X(:,5); z_dot = X(:,6);
phi = X(:,7); theta = X(:,8);
phi_dot = X(:,10); theta_dot = X(:,11); psi_dot = X(:,12);
Nrc = [0;0;-p.R];
vco = zeros(3,len);
Nomega = zeros(3,len);
vconorm = zeros(len,1);
for i = 1:len
    omega_x = theta_dot(i)*cos(phi(i))+psi_dot(i)*sin(theta(i))*sin(phi(i));
    omega_y = theta_dot(i)*sin(phi(i))-psi_dot(i)*sin(theta(i))*cos(phi(i));
    omega_z = phi_dot(i)+psi_dot(i)*cos(theta(i));
    Nomega(:,i) = [omega_x;omega_y;omega_z];
    v = [x_dot(i);y_dot(i);z_dot(i)];
    vcg = cross(Nomega(:,i),Nrc);
    vgo = v;
    vco(:,i) = vcg+vgo;
    vconorm(i) = norm(vco(:,i));
end

% ANIMATE
mode = 1; % 1 = animation, 2 = animation + record & save it
hgt_animate(X,step,p,mode)

function Xdot = bowlingballEL(t,X,p,fsymbolic)
% X = [x y z x_dot y_dot z_dot
%      phi theta psi phi_dot theta_dot psi_dot]'
% Xdot = [x_dot y_dot z_dot x_ddot y_ddot z_ddot
%         phi_dot theta_dot psi_dot phi_ddot theta_ddot psi_ddot]'

% PARAMS
g = p.g; m = p.m; R = p.R;
I1 = p.I1; I2 = p.I2; I3 = p.I3;

% COORDS
x = X(1);

% DETERMINE FRICTION
length = 18.288;
if x < .5*length
    mu = .04;
elseif x >= .5*length && x <= .7*length
    % mu = .8/length*x-.36;
    mu = 1.3/length*x-.61;
else
    mu = .3;
end

Xdot = fsymbolic(t,X,g,m,R,mu,I1,I2,I3);

end