clear all;clc;clf

N = 50;
%constants
dZ = 1/2;
dT = 2*pi/N;
A = 1/2;
T1 = .1;
H = 10; 

K = A*dT/dZ^2; %Constant alpha*dT/dZ^2

t = 0:dT:4*pi;
z = 0:dZ:3;
size_T = size(t'); size_T = size_T(1,1);
N = size(z'); N = N(1,1);

%initial conditions
for j = 1:size_T
   T(1,j) = 1 - T1*sin(t(j));
end
for i = 2:N
    T(i,1) = 1;
end

% %finite difference
% for j = 1:(size_T-1)
%     for i = 2:N
%        if(i == N) %backward diff., 1st order
%            T(i,j+1) = T(i,j) + K*(T(i,j) - 2*T(i-1,j) + T(i-2,j));
%        else %central diff, 2nd order
%            T(i,j+1) = T(i,j) + K*(T(i-1,j) - 2*T(i,j) + T(i+1,j));
%        end
%     end
%     for i = 2:N
%         if(i == N)
%             T(i,j) = T(i,j+1) + K*(T(i,j+1) - 2*T(i-1,j+1) + T(i-2,j+1));
%         else
%             T(i,j) = T(i,j+1) - K*(T(i-1,j+1) - 2*T(i,j+1) + T(i+1,j+1));
%         end
%     end
% end

%plot i = 1,2,3,5,7, z = 0,.5,1,2,3
% hold on;
% plot(t,T(1,:),'*black','LineWidth',1.5);
% plot(t,T(2,:),'xblack','LineWidth',1.5);
% plot(t,T(3,:),'+black','LineWidth',1.5);
% plot(t,T(5,:),'oblack','LineWidth',1.5);
% plot(t,T(7,:),'.black','LineWidth',1.5);
% %legend('z* = 0','z* = .5','z* = 1.0','z* = 2.0', 'z* = 3.0');
% xlabel('t*'),ylabel('T(z*,t*)');

