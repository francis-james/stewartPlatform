%% Script to linearize the symbolic equations for a Stewart platform
% Dynamic equations obtained from Lin, Lih‐Chang, and Ming‐Uei Tsay. 
% "Modeling and control of micropositioning systems using Stewart 
% platforms." Journal of Robotic Systems 17.1 (2000): 17-52.

% state = [tx ty tz thetax thetay thetaz]


%% Set up constant parameters
% Note that geometry is fixed and is not left in terms of variables
syms m Ixx Iyy Izz Ixy Ixz Iyx Iyz Izx Izy Izz;
r=1;
ang1=10*pi/180; %10 degree spacing between adjacent V spherical joint positions
ang2=110*pi/180; %110 degree spacing between 'non adjacent' spherical joint positions
baseOffsetAngle=-pi/3; % -60 degree offset for first link/universal joint on the base
g=9.806; %acceleration due to gravity
%% Set up variables
syms tx ty tz thetax thetay thetaz txd tyd tzd wx wy wz;
syms f1 f2 f3 f4 f5 f6;
syms fex fey fez tauex tauey tauez

%Rotation matrix
R=[cos(thetaz)*cos(thetay) cos(thetaz)*sin(thetay)*sin(thetax)-sin(thetaz)*cos(thetax)  ...
    cos(thetax)*sin(thetay)*cos(thetaz)+sin(thetax)*sin(thetaz); ...
    sin(thetaz)*cos(thetay) sin(thetax)*sin(thetay)*sin(thetaz)+cos(thetax)*cos(thetaz) ...
    sin(thetaz)*sin(thetay)*cos(thetax)-cos(thetaz)*sin(thetax); ...
    -sin(thetay) cos(thetay)*sin(thetax) cos(thetay)*cos(thetax)];
%Position of platform com
T=[tx;ty;tz];
%positions of leg attachments to platform
pP1=[r*cos(0); r*sin(0);0];
pP2=[r*cos(ang1); r*sin(ang1);0];
pP3=[r*cos(ang1+ang2); r*sin(ang1+ang2);0];
pP4=[r*cos(2*ang1+ang2); r*sin(2*ang1+ang2);0];
pP5=[r*cos(2*ang1+2*ang2); r*sin(2*ang1+2*ang2);0];
pP6=[r*cos(3*ang1+2*ang2); r*sin(3*ang1+2*ang2);0];
pPall=[pP1, pP2, pP3, pP4, pP5, pP6 ,pP1];
pPallX=pPall(1,:);
pPallY=pPall(2,:);
pPallZ=pPall(3,:);
plot3(pPallX,pPallY,pPallZ)

%positions of leg attachments to base
B1=[r*cos(0+baseOffsetAngle); r*sin(baseOffsetAngle); -1];
B2=[r*cos(0+baseOffsetAngle+ang2); r*sin(baseOffsetAngle+ang2); -1];
B3=[r*cos(0+baseOffsetAngle+ang2+ang1); r*sin(baseOffsetAngle+ang2+ang1); -1];
B4=[r*cos(0+baseOffsetAngle+2*ang2+ang1); r*sin(baseOffsetAngle+2*ang2+ang1); -1];
B5=[r*cos(0+baseOffsetAngle+2*ang2+2*ang1); r*sin(baseOffsetAngle+2*ang2+2*ang1); -1];
B6=[r*cos(0+baseOffsetAngle+3*ang2+2*ang1); r*sin(baseOffsetAngle+3*ang2+2*ang1); -1];
Ball=[B1, B2, B3,B4,B5,B6,B1];
BallX=Ball(1,:);
BallY=Ball(2,:);
BallZ=Ball(3,:);
hold on;
plot3(BallX,BallY,BallZ);

for i=1:6
    X=[pPall(1,i),Ball(1,i)];
    Y=[pPall(2,i),Ball(2,i)];
    Z=[pPall(3,i),Ball(3,i)];
    plot3(X,Y,Z);
    hold on;
end

for i=1:6
    BP(:,i)=T+R*pPall(:,i);
    S(:,i)=BP(:,i)-Ball(:,i);
    s(:,i)=S(:,i)/norm(S(:,i));
    Mom(:,i)=cross(R*pPall(:,i),s(:,i)); %moments
end

U=[s;Mom];
I=[Ixx -Ixy -Ixz; -Ixy Iyy -Iyz; -Ixz -Iyz Izz];
Iprime=R*I*R.';
M=[[m 0 0; 0 m 0; 0 0 m],zeros(3,3); zeros(3,3),Iprime];
% %Add expressions for c's
c45=Iprime(3,1)*wx+Iprime(3,2)*wy+Iprime(3,3)*wz;
c46=-(Iprime(2,1)*wx+Iprime(2,2)*wy+Iprime(2,3)*wz);
c54=c46;
c56=Iprime(1,1)*wx+Iprime(1,2)*wy+Iprime(1,3)*wz;
c64=-c46;
c65=-c56;
C=[zeros(3,6); zeros(3,3), [0, c45 c46; c54 0 c56; c64 c65 0]];
F=[f1;f2;f3;f4;f5;f6];
N=[fex; fey; fez+g; tauex; tauey; tauez];

%% Equation for accelerations [tx;ty;tx;thetax; thetay; thetaz]''
eq=inv(M)*(U*F-C*[txd;tyd;tzd;wx;wy;wz]);

%% Linearization
% For Jacobian
linEq=jacobian(eq,[tx ty tz thetax thetay thetaz txd tyd tzd wx wy wz F.']);

%Find initial equilibrium forces
tx=0; ty=0; tz=0; thetax=0; thetay=0; thetaz=0; txd=0; tyd=0; tzd=0; wx=0; wy=0; wz=0; Ixx=1; Iyy=1; Izz=1;
Ixy=0; Ixz=0; Iyx=0; Iyz=0; Izx=0; Izy=0;
U0=subs(U);
F0=inv(U0)*[0;0;m*g;0;0;0];

f1=F0(1);f2=F0(2);f3=F0(3);f4=F0(4);f5=F0(5);f6=F0(6);
%Evaluate at initial point
linEq0=subs(linEq);
mat=vpa(linEq0,5);

%% With Ixx, Iyy,Izz=1, other inertias=0
% [       7.8645,    5.6843e-14,  -3.4106e-13,     -0.32785,       0.97073,   2.2737e-13, 0, 0, 0, 0, 0, 0, 0.35355/m,  0.28229/m, -0.70711/m, 0.28229/m,  0.35355/m,    -0.56459/m]
% [   5.6843e-14,        7.8645,   1.1369e-13,     -0.97073,      -0.32785,  -3.4106e-13, 0, 0, 0, 0, 0, 0, 0.61237/m, -0.48895/m,          0, 0.48895/m, -0.61237/m, -9.1635e-17/m]
% [  -3.4106e-13,    1.1369e-13,       3.8829,  -2.2737e-13,    1.1369e-13,      0.65571, 0, 0, 0, 0, 0, 0, 0.70711/m,  0.82537/m,  0.70711/m, 0.82537/m,  0.70711/m,     0.82537/m]
% [   -0.32785*m,    -0.97073*m,            0,    0.23866*m, -2.2737e-13*m, 1.1369e-13*m, 0, 0, 0, 0, 0, 0,         0,    0.14332,    0.61237,   0.63227,   -0.61237,       -0.7756]
% [    0.97073*m,    -0.32785*m, 1.1369e-13*m, 4.5475e-13*m,     0.23866*m,            0, 0, 0, 0, 0, 0, 0,  -0.70711,   -0.81283,    0.35355,   0.53054,    0.35355,       0.28229]
% [ 2.2737e-13*m, -3.4106e-13*m,    0.65571*m, 2.2737e-13*m, -5.6843e-14*m,     3.2496*m, 0, 0, 0, 0, 0, 0,   0.61237,   -0.53054,    0.61237,  -0.53054,    0.61237,      -0.53054]
%  