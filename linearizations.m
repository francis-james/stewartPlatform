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
%% Set up variables
syms tx ty tz thetax thetay thetaz;

R=[cos(thetaz)*cos(thetay) cos(thetaz)*sin(thetay)*sin(thetax)-sin(thetaz)*cos(thetax)  ...
    cos(thetax)*sin(thetay)*cos(thetaz)+sin(thetax)*sin(thetaz); ...
    sin(thetaz)*cos(thetay) sin(thetax)*sin(thetay)*sin(thetaz)+cos(thetax)*cos(thetaz) ...
    sin(thetaz)*sin(thetay)*cos(thetax)-cos(thetaz)*sin(thetax); ...
    -sin(thetay) cos(thetay)*sin(thetax) cos(thetay)*cos(thetax)];

T=[tx;ty;tz];

pP1=[r*cos(0); r*sin(0);0];
pP2=[r*cos(ang1); r*sin(ang1);0];
pP3=[r*cos(ang1+ang2); r*sin(ang1+ang2);0];
pP4=[r*cos(2*ang1+ang2); r*sin(2*ang1+ang2);0];
pP5=[r*cos(2*ang1+2*ang2); r*sin(2*ang1+2*ang2);0];
pP6=[r*cos(3*ang1+2*ang2); r*sin(3*ang1+2*ang2);0];
Pall=[pP1, pP2, pP3, pP4, pP5, pP6 ,pP1];
PallX=Pall(1,:);
PallY=Pall(2,:);
PallZ=Pall(3,:);
plot3(PallX,PallY,PallZ)

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
plot3(BallX,BallY,BallZ)

for i=1:6
    X=[Pall(1,i),Ball(1,i)];
    Y=[Pall(2,i),Ball(2,i)];
    Z=[Pall(3,i),Ball(3,i)];
    plot3(X,Y,Z);
    hold on;
end

for i=1:6
    Bp(i,:)=T+R*Pall(:,i);
end

