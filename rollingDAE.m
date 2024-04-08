function fsymbolic = rollingDAE()

% PARAMS
syms t m g R mu I1 I2 I3 'real'
BI = [I1 0  0;
      0  I2 0;
      0  0  I3];

% COORDS
syms x y z x_dot y_dot z_dot x_ddot y_ddot z_ddot 'real'
syms phi theta ppsi phi_dot theta_dot psi_dot phi_ddot theta_ddot psi_ddot 'real'
n1 = [1 0 0]'; b1 = n1;
n2 = [0 1 0]'; b2 = n2;
n3 = [0 0 1]'; b3 = n3;
a = x_ddot*n1+y_ddot*n2+z_ddot*n3;

% EULER ANGLE RATES --> BODY AXIS RATES
omega_1 = phi_dot*sin(theta)*sin(ppsi)+theta_dot*cos(ppsi);
omega_2 = phi_dot*sin(theta)*cos(ppsi)-theta_dot*sin(ppsi);
omega_3 = phi_dot*cos(theta)+psi_dot;
Bomega = omega_1*b1+omega_2*b2+omega_3*b3;

% d/dt(BODY AXIS RATES) WRT BODY
alpha_1 = cos(ppsi)*theta_ddot-sin(ppsi)*theta_dot*psi_dot+sin(ppsi) ...
    *sin(theta)*phi_ddot+cos(ppsi)*sin(theta)*phi_dot*psi_dot ...
    +cos(theta)*sin(ppsi)*phi_dot*theta_dot;
alpha_2 = cos(ppsi)*sin(theta)*phi_ddot-cos(ppsi)*theta_dot*psi_dot ...
    -sin(ppsi)*theta_ddot+cos(ppsi)*cos(theta)*phi_dot*theta_dot ...
    -sin(ppsi)*sin(theta)*phi_dot*psi_dot;
alpha_3 = cos(theta)*phi_ddot+psi_ddot-sin(theta)*phi_dot*theta_dot;
Balpha = alpha_1*b1+alpha_2*b2+alpha_3*b3;

% REACTIONS
syms F_N fx fy 'real'

% COORDINATE SYSTEM TRANSFORMATION N <--> B
cphi = cos(phi); sphi = sin(phi);
cth = cos(theta); sth = sin(theta);
csi = cos(ppsi); ssi = sin(ppsi);

BRN = [cphi*csi-cth*sphi*ssi  csi*sphi+cphi*cth*ssi  ssi*sth;
      -cphi*ssi-csi*cth*sphi  cphi*csi*cth-sphi*ssi  csi*sth;
       sphi*sth               -cphi*sth                  cth];

% LMB
LMB = F_N*n3-m*g*n3 == m*a;
eqn1 = LMB(1);
eqn2 = LMB(2);
eqn3 = LMB(3);

% AMB
NM = [-R*fy;R*fx;0]; % moments in N coordinates
BM = BRN*NM; % moments in B coordinates
BHdot = BI*Balpha+cross(Bomega,BI*Bomega); % rate of change of angular momentum in B coordinates
AMB = BM == BHdot; % angular momentum balance done in B coordinates
eqn4 = AMB(1);
eqn5 = AMB(2);
eqn6 = AMB(3);

% CONSTRAINT EQUATIONS
con1 = z_ddot == 0; % normal force
con2 = x_ddot+R*sin(phi)*theta_ddot-R*cos(phi)*sin(theta)*psi_ddot ...
    +R*cos(phi)*theta_dot*phi_dot+R*sin(phi)*sin(theta)*psi_dot*phi_dot ...
    -R*cos(phi)*cos(theta)*psi_dot*theta_dot == 0; % no slip
con3 = y_ddot-R*cos(phi)*theta_ddot-R*sin(phi)*sin(theta)*psi_ddot ...
    +R*sin(phi)*theta_dot*phi_dot-R*cos(phi)*sin(theta)*psi_dot*phi_dot ...
    -R*cos(theta)*sin(phi)*psi_dot*theta_dot == 0; % no slip

% DAE MATRICIES
eqns = [eqn1 eqn2 eqn3 eqn4 eqn5 eqn6 con1 con2 con3];
vars = [x_ddot y_ddot z_ddot phi_ddot theta_ddot psi_ddot F_N fx fy];
[A,B] = equationsToMatrix(eqns,vars);

% SYMBOLIC SOL
u = A\B;
X_dot = [x_dot;y_dot;z_dot;u(1:3);phi_dot;theta_dot;psi_dot;u(4:6)];
X_sym = [x y z x_dot y_dot z_dot phi theta ppsi phi_dot theta_dot psi_dot]';
fsymbolic = matlabFunction(X_dot,'vars',{t,X_sym,g,m,R,mu,I1,I2,I3});

end

