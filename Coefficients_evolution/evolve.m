xmax = 8;
ymax = 8;
Nx = 2^8;
Ny = 2^8;
dx = 2*xmax/Nx;
dy = 2*ymax/Ny;
x = -xmax:dx:xmax-dx;
y = -ymax:dy:ymax-dy;

[X,Y] = meshgrid(x,y);

load basis.mat;

[M,~,~] = size(basis);

% evolving steps
NN = 100000;

% produce time-independent random potential
% mu = 10;
% x0 = 2;
% y0 = 2;
% bx = 2*pi*sqrt(5);
% by = 2*pi*sqrt(5);
% gr = mu*cos(bx*X/x0).*cos(by*Y/y0);

% disorder width
W = 50;
gr = (2*rand(Nx,Ny)-1)*W;

B = zeros(M,M);
for indx = 1:M
    for indy = 1:M
        B(indx,indy) = sum(sum(conj(squeeze(basis(indx,:,:))).*gr.*squeeze(basis(indy,:,:))))*dx*dy;
    end
end


% initial coefficients
u = zeros(M,NN);
u(1,1) = sqrt(1);
%u(2,1) = sqrt(0.3);
%u(3,1) = sqrt(0.2);


% eigen-values matrix
A = diag([1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9]);

% time step
h = 0.01;

for ind = 1:NN-1
    % Runge–Kutta–Fehlberg method, also called RK-45
    utemp = u(:,ind);
    k1 = h*myfunc(A,B,utemp);
    k2 = h*myfunc(A,B,utemp+k1/4);
    k3 = h*myfunc(A,B,utemp+3*k1/32+9*k2/32);
    k4 = h*myfunc(A,B,utemp+1932*k1/2197-7200*k2/2197+7296*k3/2197);
    k5 = h*myfunc(A,B,utemp+439*k1/216-8*k2+3680*k3/513-845*k4/4104);
    k6 = h*myfunc(A,B,utemp-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40);
    u(:,ind+1) = utemp + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
end
% check normalization
Norm = sum(abs(u).^2);

% reconstruct the wave function
psi0 = 0;
for ind = 1:M
    psi0 = psi0 + u(ind,1)*squeeze(basis(ind,:,:));
end

psi = 0;
for ind = 1:M
    psi = psi + u(ind,NN)*squeeze(basis(ind,:,:));
end

[Fx,Fy] = gradient(psi,dx,dy);
[Fxx,~] = gradient(Fx,dx,dy);
[~,Fyy] = gradient(Fy,dx,dy);

% compute the angular momentum, should be l, if correct
px = -1i*Fx;
py = -1i*Fy;
Lz = X.*py - Y.*px;
disp(['The numerical angular momentum is ',num2str(real(sum(sum(conj(psi).*Lz))*dx*dy))])
disp(['The exact angular momentum is ',num2str(-8)])

subplot(1,2,1)
imagesc(abs(psi0).^2);axis equal
subplot(1,2,2)
imagesc(abs(psi).^2);axis equal