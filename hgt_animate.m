function hgt_animate(X,step,p,mode)
% X = [x y z x_dot y_dot z_dot
%      phi theta psi phi_dot theta_dot psi_dot]'

% PARAMS
R = p.R;

% COORDS
x = X(:,1); y = X(:,2); z = X(:,3);
phi = X(:,7); theta = X(:,8); psi = X(:,9);

% EULER ANGLES --> DCM --> AXIS/ANGLE
len = length(phi);
ang = zeros(1,len);
a = zeros(3,len);
for i = 1:len
    cphi = cos(phi(i)); sphi = sin(phi(i));
    cth = cos(theta(i)); sth = sin(theta(i));
    csi = cos(psi(i)); ssi = sin(psi(i));

    BRN = [cphi*csi-cth*sphi*ssi  csi*sphi+cphi*cth*ssi  ssi*sth;
          -cphi*ssi-csi*cth*sphi  cphi*csi*cth-sphi*ssi  csi*sth;
           sphi*sth               -cphi*sth                  cth];

    ang(i) = acos((trace(BRN)-1)/2);
    a1 = BRN(2,3)-BRN(3,2);
    a2 = BRN(3,1)-BRN(1,3);
    a3 = BRN(1,2)-BRN(2,1);
    a(:,i) = (1/(2*sin(ang(i))).*[a1;a2;a3]);
end

% HG TRANSFORM SETUP
Q = zeros(4,4,len);
for i = 1:len
    % rotation
    if ang(i) > .01
        Qr = makehgtform('axisrotate',[a(1,i),a(2,i),a(3,i)],ang(i));
    else
        Qr = makehgtform('axisrotate',[1,0,0],0);
    end
    % translation
    Qt = makehgtform('translate',[x(i) y(i) z(i)]);
    Q(:,:,i) = Qt*Qr;
end
ax = axes('Xlim',[min(x)-R max(x)+R],'Ylim',[-0.5334-R 0.5334+R],'Zlim',[0 1]);
% ax = axes('Xlim',[min(x)-R max(x)+R],'Ylim',[min(y) max(y)],'Zlim',[0 1]);
ax.DataAspectRatio = [1 1 1];
grid on;
view(-86,10);
[xs,ys,zs] = sphere(10);
xs = xs.*R; ys = ys.*R; zs = zs.*R;
S = surface(xs,ys,zs,'FaceColor','g');
ball = hgtransform('Parent',ax);
set(S,'Parent',ball)
ball.Matrix = Q(:,:,1);
hold on
title('Bowling Ball Animation')
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

pause

% ANIMATE WITH HGTRANSFORM
if mode == 2
    % video mode
    obj = VideoWriter('bowling_ball','MPEG-4');
    obj.Quality = 75;
    obj.FrameRate = 1/step;
    open(obj);
    for i = 1:len
         ball.Matrix = Q(:,:,i);
         f = getframe(gcf);
         writeVideo(obj,f);
         drawnow
    end
    obj.close();
    pause(step)
else
    % animation mode
    for i = 1:len
         ball.Matrix = Q(:,:,i);
         drawnow
         pause(step)
    end
end

% PLOT TRAJECTORY
plot3(x,y,z)
hold off