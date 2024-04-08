function fsymbolic = slippingEL()

% PARAMS
syms t m g R mu I1 I2 I3 'real'
Nrc = [0;0;-R];

% COORDS
syms x y z x_dot y_dot z_dot x_ddot y_ddot z_ddot 'real'
syms phi theta ppsi phi_dot theta_dot psi_dot phi_ddot theta_ddot psi_ddot 'real'

% LAGRANGE MULTIPLIERS
syms lambda 'real'

% GENERALIZED COORDINATES
q = [x y z phi theta ppsi]';
q_dot = [x_dot y_dot z_dot phi_dot theta_dot psi_dot]';
q_ddot = [x_ddot y_ddot z_ddot phi_ddot theta_ddot psi_ddot]';

% LAGRANGIAN
omega_1 = phi_dot*sin(theta)*sin(ppsi)+theta_dot*cos(ppsi);
omega_2 = phi_dot*sin(theta)*cos(ppsi)-theta_dot*sin(ppsi);
omega_3 = phi_dot*cos(theta)+psi_dot;

Ek = 0.5*m*(x_dot^2+y_dot^2+z_dot^2)+0.5*(I1*omega_1^2+I2*omega_2^2+I3*omega_3^2);
Ep = m*g*z;
L = Ek - Ep;

% CONSTRAINT EQUATION
con1 = z_ddot == 0;

% EULER ANGLE RATES --> INERTIAL AXIS RATES
omega_x = theta_dot*cos(phi)+psi_dot*sin(theta)*sin(phi);
omega_y = theta_dot*sin(phi)-psi_dot*sin(theta)*cos(phi);
omega_z = phi_dot+psi_dot*cos(theta);
Nomega = [omega_x;omega_y;omega_z];

% FRICTION FORCE
vcg = cross(Nomega,Nrc);
vgo = [x_dot;y_dot;z_dot];
vco = vcg+vgo; % velocity of contact point relative to origin
voc1 = -vco./norm(vco); % direction of friction force

% GENERALIZED FORCES
% normal reaction
hz = 1;
Q_Nz = lambda*hz;
Q_N = [0 0 Q_Nz 0 0 0]';

% Rayleigh's dissipation
Q_R = -mu*q_dot;

% External force
Q_EXx = m*g*mu*voc1'*[1;0;0];
Q_EXy = m*g*mu*voc1'*[0;1;0];
Q_EX = [Q_EXx Q_EXy 0 0 0 0]';

% EOM
J = @(f,v) jacobian(f,v);
EoM = J(J(L,q_dot),q)*q_dot + J(J(L,q_dot),q_dot)*q_ddot - J(L,q)' == Q_N+Q_EX;

% SYMBOLIC SOL
eqns = [EoM;con1];
vars = [q_ddot;lambda];
[A,B] = equationsToMatrix(eqns,vars);
u = A\B;
X_dot = [x_dot;y_dot;z_dot;u(1:3);phi_dot;theta_dot;psi_dot;u(4:6)];
X_sym = [x y z x_dot y_dot z_dot phi theta ppsi phi_dot theta_dot psi_dot]';
fsymbolic = matlabFunction(X_dot,'vars',{t,X_sym,g,m,R,mu,I1,I2,I3});

end

