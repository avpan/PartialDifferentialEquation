clear all;clc;clf

N = 50;
T = 5;
X = 5;

%constants
dX = X/N;
dT = T/N;
v0 = 0; %velocity
rho_0 = 5; %density
p0 = 1; %gas pressure
gamma = 1.2; %adiabadtic index
cs0 = gamma*p0/rho_0;

%perturbations
prime1 = cs0*1E-10;
prime2 = rho_0*1E-10;
prime3 = p0*1E-10;

K1 = dT/dX/rho_0; %Constants 
K2 = rho_0*dT/dX;
K = dT/dX;
t = 0:dT:T;
x = 0:dX:X;
size_T = length(t);
N = length(x);

%initial conditions
for j = 1:size_T
   v(1,j) = v0 + prime1*j;
   rho(1,j) = rho_0 + prime2*j;
   p(1,j) = p0 + prime3*j;
   
   v(N,j) = v0;
   rho(N,j) = rho_0;
   p(N,j) = p0;
end

for i = 1:N
    v(i,1) = v0 + prime1*i;
    rho(i,1) = rho_0 + prime2*i;
    p(i,1) = p0 + prime3*i;
end

%finite difference
for j = 1:(size_T-1)
    for i = 2:N-1
        v(i,j+1) = v(i,j) - K1*(p(i+1,j) - p(i,j));
        rho(i,j+1) = rho(i,j) - K2*(v(i+1,j) - v(i,j));
        p(i,j+1) = cs0^2*rho(i,j+1);
    end
end

% 
% %finite difference
% for j = 1:(size_T-1)
%     for i = 2:N-1
%            v(i,j+1) = v(i,j) - K*(v(i+1,j) - v(i,j)) - K*(p(i+1,j)-p(i,j));
%            rho(i,j+1) = rho(i,j) - K*(rho(i+1,j)*v(i+1,j) - rho(i,j)*v(i,j));
%            p(i,j+1) = p0*(rho(i,j+1)/rho_0)^gamma;
%     end
% end

colormap('white');
figure(1);
surf(x,t,v);
set(gca, 'XDir', 'reverse');
xlabel('x distance');ylabel('Time (t)');zlabel('Fluid Velocity (v)');

colormap('white');
figure(2);
surf(rho,t,v);
xlabel('gas density (\rho)');ylabel('Time (t)');zlabel('Fluid Velocity (v)');

colormap('white');
figure(3);
surf(p,t,v);
xlabel('gas pressure (p)');ylabel('Time (t)');zlabel('Fluid Velocity (v)');

colormap('white');
% figure(4);
% surf(p,t,rho);
% xlabel('gas pressure (p)');ylabel('Time (t)');zlabel('gas density (\rho)');
% colormap('white');
