function fsymbolic = rollingEL()

% PARAMS
syms t m g R mu I1 I2 I3 'real'

% COORDS
syms x y z x_dot y_dot z_dot x_ddot y_ddot z_ddot 'real'
syms phi theta ppsi phi_dot theta_dot psi_dot phi_ddot theta_ddot psi_ddot 'real'

% LAGRANGE MULTIPLIERS
syms lambda_f lambda_g lambda 'real'

% LAGRANGIAN
omega_1 = phi_dot*sin(theta)*sin(ppsi)+theta_dot*cos(ppsi);
omega_2 = phi_dot*sin(theta)*cos(ppsi)-theta_dot*sin(ppsi);
omega_3 = phi_dot*cos(theta)+psi_dot;

Ek = 0.5*m*(x_dot^2+y_dot^2+z_dot^2)+0.5*(I1*omega_1^2+I2*omega_2^2+I3*omega_3^2);
Ep = m*g*z;
L = Ek - Ep;

% CONSTRAINT EQUATIONS (NO SLIP)
con1 = x_ddot+R*sin(phi)*theta_ddot-R*cos(phi)*sin(theta)*psi_ddot ...
    +R*cos(phi)*theta_dot*phi_dot+R*sin(phi)*sin(theta)*psi_dot*phi_dot ...
    -R*cos(phi)*cos(theta)*psi_dot*theta_dot == 0;
con2 = y_ddot-R*cos(phi)*theta_ddot-R*sin(phi)*sin(theta)*psi_ddot ...
    +R*sin(phi)*theta_dot*phi_dot-R*cos(phi)*sin(theta)*psi_dot*phi_dot ...
    -R*cos(theta)*sin(phi)*psi_dot*theta_dot == 0;
con3 = z_ddot == 0;

% GENERALIZED FORCES
fx = 1;
fy = 0;
fz = 0;
fphi = 0;
ftheta = R*sin(phi);
fpsi = -R*sin(theta)*cos(phi);

gx = 0;
gy = 1;
gz = 0;
gphi = 0;
gtheta = -R*cos(phi);
gpsi = -R*sin(theta)*sin(phi);

hz = 1;

Qx = lambda_f*fx+lambda_g*gx;
Qy = lambda_f*fy+lambda_g*gy;
Qz = lambda_f*fz+lambda_g*gz+lambda*hz;
Qphi = lambda_f*fphi+lambda_g*gphi;
Qtheta = lambda_f*ftheta+lambda_g*gtheta;
Qpsi = lambda_f*fpsi+lambda_g*gpsi;

Q = [Qx Qy Qz Qphi Qtheta Qpsi]';

% GENERALIZED COORDS
q = [x y z phi theta ppsi]';
q_dot = [x_dot y_dot z_dot phi_dot theta_dot psi_dot]';
q_ddot = [x_ddot y_ddot z_ddot phi_ddot theta_ddot psi_ddot]';

% EOM
J = @(f,x) jacobian(f,x);
EoM = J(J(L,q_dot),q)*q_dot + J(J(L,q_dot),q_dot)*q_ddot - J(L,q)' == Q;

% SYMBOLIC SOL
eqns = [EoM;con1;con2;con3];
vars = [q_ddot;lambda_f;lambda_g;lambda];
[A,B] = equationsToMatrix(eqns,vars);
u = A\B;
X_dot = [x_dot;y_dot;z_dot;u(1:3);phi_dot;theta_dot;psi_dot;u(4:6)];
X_sym = [x y z x_dot y_dot z_dot phi theta ppsi phi_dot theta_dot psi_dot]';
fsymbolic = matlabFunction(X_dot,'vars',{t,X_sym,g,m,R,mu,I1,I2,I3});
end

