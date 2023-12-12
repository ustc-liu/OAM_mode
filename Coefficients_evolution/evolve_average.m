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
% disorder width
W = 100;
Nsamples = 10;

% eigen-values matrix
A = diag([1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9]);
% time step
h = 0.01;
psi2 = zeros(M,NN);
for indn = 1:Nsamples
    indn
    dV = (2*rand(Nx,Ny)-1)*W;
    B = zeros(M,M);
    for indx = 1:M
        for indy = 1:M
            B(indx,indy) = sum(sum(conj(squeeze(basis(indx,:,:))).*dV.*squeeze(basis(indy,:,:))))*dx*dy;
        end
    end

    % initial coefficients
    u = zeros(M,NN);
    u(1,1) = sqrt(10)*exp(1i*rand*2*pi);
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
    psi2 = psi2 + abs(u).^2;
end
psi2 = psi2/Nsamples;