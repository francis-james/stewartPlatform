%% Dynamic simulation of tewart platform
%%state=[tx ty tz thetax thetay thetaz txd tyd tzd wx wy wz]
state0=zeros(12,1);
state0(3)=1.1;
stewart=StewartPlatform(state0);
U0=stewart.get_U(stewart.x);
F0=inv(U0)*[0;0;stewart.m*stewart.g;0;0;0];
[T,Y]=ode45(@(t,y)dynamicsPlatform(t,y,stewart,F0), [0 20], stewart.x);


fig=figure();
axs=axes('Parent',fig);
for i = 1:10:length(T)
    stewart.plot(Y(i,:)',fig,axs)
    hold off
    pause(0.05)
end
