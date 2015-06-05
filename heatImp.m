clear all;clc;clf;

N = 64;
%constants
dZ = 1/2;
dT = 2*pi/N;
A = 1/2;
T1 = .1;
H = 10; 

s = A*dT/dZ^2; %Constant alpha*dT/dZ^2

t = 0:dT:4*pi;
z = 0:dZ:10;
size_T = length(t);
N = length(z);

A = sparse(N,N);
for i = 2:N-1
    A(i,i-1) = -s;
    A(i,i) = 1+2*s;
    A(i,i+1) = -s;
end

%set boundary conditions
A(1,1) = 1;
A(N,N) = 1;

B = zeros(N,1);
B(:,1) = 1;
B(1) = 1; B(N) = 1;
B;
T(:,1) = A\B;
B = A\B;

%implicit solve for AT = B
for i = 2:size_T
    B(1) = 1+T1*sin(t(i)); 
    B(N) = 1;
    
    T(:,i) = A\B;
    B = A\B;
end

%plot i = 1,2,3,5,7, z = 0,.5,1,2,3
hold on;
plot(t,T(1,:),'*black','LineWidth',1.5);
plot(t,T(2,:),'xblack','LineWidth',1.5);
plot(t,T(3,:),'+black','LineWidth',1.5);
plot(t,T(5,:),'oblack','LineWidth',1.5);
plot(t,T(7,:),'.black','LineWidth',1.5);
legend('z* = 0','z* = .5','z* = 1.0','z* = 2.0', 'z* = 3.0');
xlabel('t*'),ylabel('T(z*,t*)');

% hold on;
% plot(t,T(1,:),'*black');
% for i = 7:7
%     plot(t,T(i,:),'oblack');
% end
% legend('z* = 0','z* = .5');
% xlabel('t*'),ylabel('T(z*,t*)');
