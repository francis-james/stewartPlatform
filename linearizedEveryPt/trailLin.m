%% Trial 1 for getting linearized function 
stewart=StewartPlatform(zeros(12,1));
state=stewart.x;
U0=stewart.get_U(state);
f=inv(U0)*[0;0;stewart.m*stewart.g;0;0;0];
linEq=getLinEq(state,stewart,f);