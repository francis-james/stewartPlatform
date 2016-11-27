clear
close all

sp = StewartPlatform(zeros(18,1));
B0 = sp.get_B(zeros(6,1));
b = B0(7:12,:);
u0 = inv(b)*[0;0;-9.806;0;0;0];  %linearize about point that cancels g
%[A,B]=sp.linear_f(zeros(18,1),u0);
% [A,B]=sp.linear_full_f(zeros(18,1),u0);
% C = sp.linear_full_h(zeros(18,1));
% 
% 
% sys = ss(A(1:12,1:12),B(1:12,:),eye(12),zeros(12,6));
% 
% Q = diag([10 10 10 10 10 10 5 5 5 5 5 5]);
% R = diag([.1 .1 .1 .1 .1 .1]);
% 
% [K,S,E] = lqr(sys,Q,R);

% blah = load('C:\Users\tapgar\Documents\MATLAB\stewartPlatform_francis\stewartPlatform-master\K.mat');
T=10;
blah=load('K.mat');
K = -blah.K;

% K =[ -0.0854    1.4559   -0.2199    0.5890    2.0053   -0.7272   -0.5343   -3.0301   -1.0041    1.0523    1.8227   -1.5622;...
%    -0.5298   -1.0909   -0.6660   -1.0067    2.8033    0.8879   -0.0000    2.6360   -1.1590   -1.3180    2.2828    1.3383;...
%    -1.2181   -0.8019   -0.2199   -2.0311   -0.4925   -0.7272    2.8913    1.0524   -1.0041   -2.1047   -0.0000   -1.5622;...
%     1.2097    0.0866   -0.6660   -1.9244   -2.2735    0.8879   -2.2828   -1.3180   -1.1590   -1.3180   -2.2828    1.3383;...
%     1.3035   -0.6540   -0.2199    1.4421   -1.5127   -0.7272   -2.3570    1.9778   -1.0042    1.0524   -1.8227   -1.5622;...
%    -0.6799    1.0043   -0.6660    2.9311   -0.5298    0.8879    2.2828   -1.3180   -1.1590    2.6360   -0.0000    1.3383];
% K = -K;
x = zeros(12,1);
% x(2,1) = -0.05;
% x(3,1) = 1.0;
% x(4,1) = 0.5;
% x(6,1) = pi/8;
% x(5,1) = 1;
fe = zeros(6,1);


dt = 0.01;
N=floor(T/dt);
t=(0:dt:T);
x_hist = zeros(N,18);
xe_hist = zeros(N,18);
xode_hist=[];
link_forces = zeros(N,6);
xm = x;

for i = 1:1:N
    fe=getFe(t(i));
    %u = zeros(6,1);
    u = -K*xm + u0;
    link_forces(i,:) = u';
    [time,xode]=ode45(@(t,y)dynamicsPlatform(t,y,u,sp),[t(i) t(i)+dt],[x;fe]);
%     xd = sp.f([x;fe],u);
%     x = x + xd(1:12,1)*dt;
    x=xode(end,1:12).';
    xode_hist=[xode_hist;xode(end,1:12),fe.'];

    z = sp.SimulateMeasurement([x;fe]);
    sp = sp.UpdateEKF(u,z,dt);
    x_hist(i,:) = [x;fe]';
    xe_hist(i,:) = sp.xest';
    xm = sp.xest(1:12,1);
    
    tau = sp.get_Torque(sp.xest);
    u0 = -inv(b)*[-sp.xest(13,1);-sp.xest(14,1);-sp.xest(15,1)-9.806;-tau];
    
end


for i = 1:10:length(x_hist)
    sp.plot(x_hist(i,:)')
    hold off
    pause(0.05)
end

i = linspace(0,10,N);

figure;
subplot(2,1,1)
plot(i,x_hist(:,1),'b-.')
hold on
plot(i,x_hist(:,2),'r-.')
plot(i,x_hist(:,3),'g-.')
plot(i,xe_hist(:,1),'b')
plot(i,xe_hist(:,2),'r')
plot(i,xe_hist(:,3),'g')
legend('x','y','z','xe','ye','ze')

subplot(2,1,2)
plot(i,x_hist(:,4),'b-.')
hold on
plot(i,x_hist(:,5),'r-.')
plot(i,x_hist(:,6),'g-.')
plot(i,xe_hist(:,4),'b')
plot(i,xe_hist(:,5),'r')
plot(i,xe_hist(:,6),'g')
legend('px','py','pz','pxe','pye','pze')

figure;
subplot(2,1,1)
plot(i,x_hist(:,13),'b-.')
hold on
plot(i,x_hist(:,14),'r-.')
plot(i,x_hist(:,15),'g-.')
plot(i,xe_hist(:,13),'b')
plot(i,xe_hist(:,14),'r')
plot(i,xe_hist(:,15),'g')
legend('fx','fy','fz','fxe','fye','fze')

subplot(2,1,2)
plot(i,x_hist(:,16),'b-.')
hold on
plot(i,x_hist(:,17),'r-.')
plot(i,x_hist(:,18),'g-.')
plot(i,xe_hist(:,16),'b')
plot(i,xe_hist(:,17),'r')
plot(i,xe_hist(:,18),'g')
legend('Tx','Ty','Tz','Txe','Tye','Tze')

figure;
plot(i,link_forces)

figure;
plot(i,x_hist(:,7:12))

hold on
plot(i,xe_hist(:,7:12))