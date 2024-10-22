clear
close all

sp = StewartPlatform(zeros(18,1));
B0 = sp.get_B(zeros(6,1));
b = B0(7:12,:);
u0 = inv(b)*[0;0;-9.806;0;0;0];  %linearize about point that cancels g

% blah = load('C:\Users\tapgar\Documents\MATLAB\stewartPlatform_francis\stewartPlatform-master\K.mat');
T=20;
Kmat=load('K.mat');
K = -Kmat.K;


x = zeros(12,1); %[tx,ty,tz,thetax,thetay,thetaz,txd,tyd,tzd,wx,wy, wz].'
fe = zeros(6,1); %[fx, fy, fz, Mx, My, Mz].'

dt = 0.01; %discretized time interval for Kalman updates
N=floor(T/dt);
t=(0:dt:T);
x_hist = zeros(N,18);
xe_hist = zeros(N,18);
% xode_hist=[];
link_forces = zeros(N,6);
xm = x;

for i = 1:1:N
    fe=getFe(t(i));
%     u = zeros(6,1);
    u = -K*xm + u0;
    link_forces(i,:) = u';
    [time,xode]=ode45(@(t,y)dynamicsPlatform(t,y,u,sp),[t(i) t(i)+dt],[x;fe]);
    x=xode(end,1:12).';
%     xode_hist=[xode_hist;xode(end,1:12),fe.'];

    z = sp.SimulateMeasurement([x;fe]);
    sp = sp.UpdateEKF(u,z,dt);
    x_hist(i,:) = [x;fe]';
    xe_hist(i,:) = sp.xest';
    xm = sp.xest(1:12,1);
    
    tau = sp.get_Torque(sp.xest);
%     u0=zeros(6,1);
    u0 = -inv(b)*[-sp.xest(13,1);-sp.xest(14,1);-sp.xest(15,1)-9.806;-tau];
    
end

% fig=figure();
% axs=axes('Parent',fig);
v=VideoWriter('fixedPtForce2.avi');
v.FrameRate=15;
open(v);
for i = 1:10:length(x_hist)
    sp.plot(x_hist(i,:)')
    hold off
    pause(0.05)
    frame=getframe(gcf);
    writeVideo(v,frame);
end
close(v);

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