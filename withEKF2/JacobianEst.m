function [ F ] = JacobianEst( func, x )

    fx=feval(func,x);
    n=length(fx);
    %F = zeros(n,n);
    eps=[-0.01,-0.005,-0.001,0.001,0.005,0.01]; % could be made better
    xperturb=x;
    m = n;
    if (length(x) == 18)
        m = length(x);
    end
    for i=1:m
        temp = zeros(n,length(eps));
        for j = 1:1:length(eps)
            xperturb(i)=xperturb(i)+eps(j);
            temp(:,j)=(feval(func,xperturb)-fx)./eps(j);
            xperturb(i)=x(i);
        end
        F(:,i) = mean(temp,2);
    end

end

