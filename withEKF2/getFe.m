function fe=getFe(t)
fe=zeros(6,1);
fe(4,1) = 1*sin(t*0.01);

fe(1:3,1) = zeros(3,1);
if (t > 3)
    fe(1,1) = 1.0+ randn(1)*0.01;
    fe(2,1) = -0.5 + randn(1)*0.01;
    fe(3,1) = -10.0+ randn(1)*0.01;
end