%% Check for controllability
m=1;
mat=[       7.8645,    5.6843e-14,  -3.4106e-13,     -0.32785,       0.97073,   2.2737e-13, 0, 0, 0, 0, 0, 0, 0.35355/m,  0.28229/m, -0.70711/m, 0.28229/m,  0.35355/m,    -0.56459/m;...
   5.6843e-14,        7.8645,   1.1369e-13,     -0.97073,      -0.32785,  -3.4106e-13, 0, 0, 0, 0, 0, 0, 0.61237/m, -0.48895/m,          0, 0.48895/m, -0.61237/m, -9.1635e-17/m;...
  -3.4106e-13,    1.1369e-13,       3.8829,  -2.2737e-13,    1.1369e-13,      0.65571, 0, 0, 0, 0, 0, 0, 0.70711/m,  0.82537/m,  0.70711/m, 0.82537/m,  0.70711/m,     0.82537/m;...
   -0.32785*m,    -0.97073*m,            0,    0.23866*m, -2.2737e-13*m, 1.1369e-13*m, 0, 0, 0, 0, 0, 0,         0,    0.14332,    0.61237,   0.63227,   -0.61237,       -0.7756;...
    0.97073*m,    -0.32785*m, 1.1369e-13*m, 4.5475e-13*m,     0.23866*m,            0, 0, 0, 0, 0, 0, 0,  -0.70711,   -0.81283,    0.35355,   0.53054,    0.35355,       0.28229;...
 2.2737e-13*m, -3.4106e-13*m,    0.65571*m, 2.2737e-13*m, -5.6843e-14*m,     3.2496*m, 0, 0, 0, 0, 0, 0,   0.61237,   -0.53054,    0.61237,  -0.53054,    0.61237,      -0.53054];
 
a1=[zeros(6,6),eye(6)];
a2= mat(:,1:12);
A=[a1;a2];

b2=mat(:,13:18);
b1=zeros(size(b2));
B=[b1;b2];
P=[B A*B];
res=rref(P);

%% Check svd to make sure it's not nearly singular
[U,S,V]=svd(P);
c=cond(P); %returns ratio of largest singular value to smallest

if c>10
    display('System is not controllable');
    exit;
end

display('Huzzah! It''s controllable');
%% From the P matrix, mu1=mu2=...=mu6=2
mu=[2;2;2;2;2;2];
M=[P(:,1),P(:,7),P(:,2),P(:,8),P(:,3),P(:,9),P(:,4),P(:,10),P(:,5),P(:,11),P(:,6),P(:,12)];
Minv=inv(M);
T=[Minv(2,:);Minv(2,:)*A;Minv(4,:);Minv(4,:)*A;Minv(6,:);Minv(6,:)*A;Minv(8,:);Minv(8,:)*A;Minv(10,:);Minv(10,:)*A;Minv(12,:);Minv(12,:)*A];

Abar=T*A*inv(T);
Bbar=T*B;

Asubdes=[0 1; -6 -5]; %placing poles at -5,-5 for each subsystem
Adesired=zeros(size(A));
startInd=1; endInd=startInd+mu(1)-1;
for i=1:length(mu)
    Adesired(startInd:endInd,startInd:endInd)=Asubdes;
    if i+1<=length(mu)
        startInd=startInd+mu(i);
        endInd=startInd+mu(i+1)-1;
    end
end
Kbar=pinv(Bbar)*(Adesired-Abar);
K=Kbar*T;