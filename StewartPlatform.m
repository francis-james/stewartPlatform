classdef StewartPlatform
    %Author: Taylor Apgar
    %STEWARTPLATFORM Summary of this class goes here
    %   Detailed explanation goes here
    properties
        
        I = eye(3);
        
        g=-9.806; %acceleration due to gravity
        m = 1;
        h = 0.1;
        
        x = zeros(12,1);
        
        r
        pP_p
        pB_b

    end
    
    methods
        
      function obj = StewartPlatform(x_init)
        obj.x = x_init;
        r = 1;
        obj.r = r;
%         obj.I(1,1) = (1/12)*obj.m*(3*obj.r^2 + obj.h^2);
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
        R=obj.get_R(obj.x(4:6));
        obj.pP_p=R*[pP1, pP2, pP3, pP4, pP5, pP6];
        
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
      
      function H = get_H(obj, state)
          q=state(1:6);
          H = zeros(6);
          H(1,1) = obj.m;
          H(2,2) = obj.m;
          H(3,3) = obj.m;
          
          R = obj.get_R(q(4:6,1));
          H(4:6,4:6) = R*obj.I*R';
      end
      
      function C = get_C(obj, state)
          q=state(1:6);
          qd=state(7:12);
          C = zeros(6);
          R = obj.get_R(q(4:6,1));
          RdRt = [0 -qd(6,1) qd(5,1); qd(6,1) 0 -qd(4,1); -qd(5,1) qd(4,1) 0];
          C(4:6, 4:6) = RdRt*R*obj.I*R';
      end
      
      function G = get_G(obj, state, fe, fp)
          q=state(1:6);
          G = zeros(6,1);
          G(1:3,1) = obj.m*[0;0;obj.g] + fe;
          R = obj.get_R(q(4:6,1));
          G(4:6,1) = cross((fp-q(1:3)),fe);
      end
      
      function U = get_U(obj, state)
          q=state(1:6);
          % B = U <= [S/|S| M]
%           U = zeros(18,6);
          R = obj.get_R(q(4:6,1));
          
          S = q(1:3,1)*ones(1,6) + R*obj.pP_p - obj.pB_b;
          for i = 1:1:6
              S(:,i) = S(:,i)./norm(S(:,i));
              M(:,i) = cross(R*obj.pP_p(:,i),S(:,i));
          end
          
          U = [S; M];
      end
      
%       function xd = f(obj,x,u)
%           xd = zeros(18,1);
%           xd(1:6,1) = x(7:12,1);
%           full_B = obj.get_B(x(1:6,1));
%           xd(7:12,1) = inv(obj.get_H(x(1:6,1)))*(full_B(7:12,:)*u...
%               - obj.get_C(x(1:6,1),x(7:12,1))*x(7:12,1) - obj.get_G(x(1:6,1),fe,fp));
%       end
%       
%       function xd = stackedf(obj,x)
%           u = x(19:end,1);
%           x = x(1:18,1);
%           xd = obj.f(x,u);
%       end
      
      function [A,B] = linear_f(obj,x,u)
          [AB, err] = jacobianest(@obj.stackedf,[x;u]);
          A = AB(:,1:18);
          B = AB(:,19:end);
      end
      
      function plot(obj,x, fig, axs)
          
          if ~exist('fig','var')
                fig=figure();
          end

          if ~exist('axs','var')
                axs=gca;
          end
            
          BallX=[obj.pB_b(1,:),obj.pB_b(1,1)];
          BallY=[obj.pB_b(2,:),obj.pB_b(2,1)];
          BallZ=[obj.pB_b(3,:),obj.pB_b(3,1)];
          
          plot3(BallX,BallY,BallZ, 'Parent',axs)
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
          plot3(pPallX,pPallY,pPallZ,'Parent',axs)
          
          xlim([-1.5 1.5])
          ylim([-1.5 1.5])
          zlim([-1 1.5])
      end
        
    end
    
end


