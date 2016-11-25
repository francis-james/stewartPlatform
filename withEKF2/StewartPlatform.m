classdef StewartPlatform
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        I = eye(3);
        
        g=9.806; %acceleration due to gravity
        m = 1;
        height = 0.1;
        
        x = zeros(18,1);
        
        r
        pP_p
        pB_b
        
        % for EKF
        model_C = 0.0001.*eye(18);
        meas_C = 0.00000001.*eye(9);
        
        P = eye(18);
        xest

    end
    
    methods
        
      function obj = StewartPlatform(x_init)
        obj.x = x_init;
        obj.xest = x_init;
        
        
        obj.model_C(13,13) = 10;
        obj.model_C(14,14) = 10;
        obj.model_C(15,15) = 10;
        obj.model_C(16,16) = 1;
        obj.model_C(17,17) = 1;
        obj.model_C(18,18) = 1;
        
        r = 1;
        obj.r = r;
%         obj.I(1,1) = (1/12)*obj.m*(3*obj.r^2 + obj.height^2);
%         obj.I(2,2) = obj.I(1,1);
%         obj.I(3,3) = 0.5*obj.m*obj.r^2;
        
        ang1=10*pi/180; %10 degree spacing between adjacent V spherical joint positions
        ang2=110*pi/180; %110 degree spacing between 'non adjacent' spherical joint positions
        baseOffsetAngle=-pi/3; % -60 degree offset for first link/universal joint on the base
        %positions of leg attachments to platform
        pP1=[r*cos(0); r*sin(0);0];
        pP2=[r*cos(ang1); r*sin(ang1);0];
        pP3=[r*cos(ang1+ang2); r*sin(ang1+ang2);0];
        pP4=[r*cos(2*ang1+ang2); r*sin(2*ang1+ang2);0];
        pP5=[r*cos(2*ang1+2*ang2); r*sin(2*ang1+2*ang2);0];
        pP6=[r*cos(3*ang1+2*ang2); r*sin(3*ang1+2*ang2);0];
        obj.pP_p=[pP1, pP2, pP3, pP4, pP5, pP6];
        
        B1=[r*cos(0+baseOffsetAngle); r*sin(baseOffsetAngle); -1];
        B2=[r*cos(0+baseOffsetAngle+ang2); r*sin(baseOffsetAngle+ang2); -1];
        B3=[r*cos(0+baseOffsetAngle+ang2+ang1); r*sin(baseOffsetAngle+ang2+ang1); -1];
        B4=[r*cos(0+baseOffsetAngle+2*ang2+ang1); r*sin(baseOffsetAngle+2*ang2+ang1); -1];
        B5=[r*cos(0+baseOffsetAngle+2*ang2+2*ang1); r*sin(baseOffsetAngle+2*ang2+2*ang1); -1];
        B6=[r*cos(0+baseOffsetAngle+3*ang2+2*ang1); r*sin(baseOffsetAngle+3*ang2+2*ang1); -1];
        obj.pB_b=[B1, B2, B3,B4,B5,B6];
               
      end
      
      function R = get_R(obj,t)
         Rx = [1 0 0; 0 cos(t(1)) -sin(t(1)); 0 sin(t(1)) cos(t(1))];
         Ry = [cos(t(2)) 0 sin(t(2)); 0 1 0; -sin(t(2)) 0 cos(t(2))];
         Rz = [cos(t(3)) -sin(t(3)) 0; sin(t(3)) cos(t(3)) 0; 0 0 1];
         R = Rx*Ry*Rz;
      end
      
      function M = get_M(obj, q)
          M = zeros(6);
          M(1,1) = obj.m;
          M(2,2) = obj.m;
          M(3,3) = obj.m;
          
          R = obj.get_R(q(4:6,1));
          M(4:6,4:6) = R*obj.I*R';
      end
      
      function C = get_C(obj, q, qd)
          C = zeros(6);
          R = obj.get_R(q(4:6,1));
          RdRt = [0 -qd(6,1) qd(5,1); qd(6,1) 0 -qd(4,1); -qd(5,1) qd(4,1) 0];
          C(4:6, 4:6) = RdRt*R*obj.I*R';
      end
      
      function G = get_G(obj, q, fe, fp)
          G = zeros(6,1);
          G(1:3,1) = obj.m*[0;0;obj.g] + fe;
          R = obj.get_R(q(4:6,1));
          G(4:6,1) = cross(R*fp,fe);
      end
      
      function tau = get_Torque(obj, x)
         R = obj.get_R(x(4:6,1));
         tau = cross(R*x(16:18,1),x(13:15,1));
      end
      
      function B = get_B(obj, q)
         
          % B = U <= [S/|S| M]
          B = zeros(18,6);
          R = obj.get_R(q(4:6,1));
          
          S = q(1:3,1)*ones(1,6) + R*obj.pP_p - obj.pB_b;
          for i = 1:1:6
              S(:,i) = S(:,i)./norm(S(:,i));
              M(:,i) = cross(R*obj.pP_p(:,i),S(:,i));
          end
          
          B(7:12,:) = [S; M];
      end
      
      function hx = h(obj,x) %measurement model
          hx = zeros(9,1);
          
          R = obj.get_R(x(4:6,1));
          S = x(1:3,1)*ones(1,6) + R*obj.pP_p - obj.pB_b;
          
          for i = 1:1:6
              hx(i,1) = norm(S(:,i));
          end
          hx(7:9,1) = x(16:18,1);
          
      end
      
      function xd = f(obj,x,u)
          xd = zeros(18,1);
          xd(1:6,1) = x(7:12,1);
          full_B = obj.get_B(x(1:6,1));
          xd(7:12,1) = inv(obj.get_M(x(1:6,1)))*(full_B(7:12,:)*u...
              - obj.get_C(x(1:6,1),x(7:12,1))*x(7:12,1) - obj.get_G(x(1:6,1),x(13:15,1),x(16:18,1)));
      end
      
      function xd = stackedf(obj,x)
          u = x(19:end,1);
          x = x(1:18,1);
          xd = obj.f(x,u);
      end
      
      function [A,B] = linear_f(obj,x,u)
          AB = JacobianEst(@obj.stackedf,[x;u]);
          A = AB(:,1:18);
          B = AB(:,19:end);
      end
      
      function xd = stackedfull(obj,x)
          u = x(19:end,1);
          x = x(1:18,1);
          xd = obj.f(x,u);
      end
      
      function [A,B] = linear_full_f(obj,x,u)
          AB = jacobianest(@obj.stackedfull,[x;u]);
          A = AB(:,1:18);
          B = AB(:,19:end);
      end
      
      function H = linear_full_h(obj,x)
         H = jacobianest(@obj.h,x); 
      end
      
      function H = linear_h(obj,x)
          H = JacobianEst(@obj.h,x);
      end
      
      function z = SimulateMeasurement(obj,x)
         
          z = obj.h(x); %add noise later
          z(1:6,1) = z(1:6,1);% + randn(6,1).*0.001;
          z(7:9,1) = z(7:9,1) + randn(3,1).*0.001;
          
      end
      
      function obj=UpdateEKF(obj, u, z, dt)
          obj=obj.ModelUpdate(u,dt);
          obj=obj.MeasUpdate(z);
      end
        
      function obj= ModelUpdate(obj, u, dt)
            xd = obj.f(obj.xest,u);
            obj.xest = obj.xest + xd.*dt;
            A = eye(18) + (dt.*obj.linear_f(obj.xest,u));
            obj.P = A*obj.P*A' + obj.model_C;
        end

        function obj = MeasUpdate(obj, z)
            hx = obj.h(obj.xest);
            H = obj.linear_h(obj.xest);
            K = (obj.P*H')/(H*obj.P*H' + obj.meas_C);
            obj.xest = obj.xest + K*(z-hx);
            obj.P = (eye(18) - K*H)*obj.P;
        end
      
      
      
      
      
      function plot(obj,x)
          BallX=[obj.pB_b(1,:),obj.pB_b(1,1)];
          BallY=[obj.pB_b(2,:),obj.pB_b(2,1)];
          BallZ=[obj.pB_b(3,:),obj.pB_b(3,1)];
          
          plot3(BallX,BallY,BallZ)
          for i=1:6
              BP(:,i) = x(1:3,1) + obj.get_R(x(4:6,1))*obj.pP_p(:,i);
          end
          hold on;
          for i=1:1:6
              X=[BP(1,i),obj.pB_b(1,i)];
              Y=[BP(2,i),obj.pB_b(2,i)];
              Z=[BP(3,i),obj.pB_b(3,i)];
              plot3(X,Y,Z);
          end
          
          pPallX=[BP(1,:),BP(1,1)];
          pPallY=[BP(2,:),BP(2,1)];
          pPallZ=[BP(3,:),BP(3,1)];
          plot3(pPallX,pPallY,pPallZ)
          
          xlim([-1.5 1.5])
          ylim([-1.5 1.5])
          zlim([-1 1.5])
      end
        
    end
    
end

