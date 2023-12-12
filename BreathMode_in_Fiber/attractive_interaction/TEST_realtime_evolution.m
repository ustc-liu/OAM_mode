xmax = 8;
ymax = 8;
Nx = 2^7;
Ny = 2^7;
dx = 2*xmax/Nx;
dy = 2*ymax/Ny;
x = -xmax:dx:xmax-dx;
y = -ymax:dy:ymax-dy;

%x = [0:dx:(xmax-dx) -xmax:dx:-dx];
%y = x;

% space-spatial meshgrid
[X,Y] = meshgrid(x,y);
[phi,r] = cart2pol(X,Y);

% momentum-space meshgrid
kx = [0:Nx/2-1 -Nx/2:-1]*pi/xmax;
ky = [0:Ny/2-1 -Ny/2:-1]*pi/ymax;

[skx,sky] = meshgrid(kx,ky);
dkx = pi/xmax;
dky = pi/ymax;

T = 0.5*skx.^2 + 0.5*sky.^2;

Vr = 0.5*r.^2;

% Suzuki-Trotter coefficients
s = 1/(2-2^(1/3));
p1 = s/2;
p2 = s;
p3 = (1-s)/2;
p4 = (1-2*s);
p5 = (1-s)/2;
p6 = s;
p7 = s/2;

Nt = 100000;

% imagnary-time evolution, time step
dt = 1i*0.001;
% initial wave-function
load GS.mat;

A = exp(-p7*T*dt);
B = exp(-p3*T*dt);

indm = 1;

delta_r = zeros(1,Nt/1000);

for ind = 0:Nt-1
    % evolving in momentum space
    Psi_k = fft2(Psi);
    TA = A.*Psi_k;
    
    % evolving in spatial space
    FA = ifft2(TA);
    % non-linear and potential part
    u = Vr - abs(FA).^2;
    Psi_x = exp(-p6*dt*u).*FA;
    
    % evolving in momentum space
    Psi_k = fft2(Psi_x);
    TA = B.*Psi_k;
    
    % evolving in momentum space
    FA = ifft2(TA);
    u = Vr - abs(FA).^2;
    Psi_x = exp(-p4*dt*u).*FA;
    
    % evolving in spatial space
    Psi_k = fft2(Psi_x);
    TA = B.*Psi_k;
    
    % evolving in momentum space
    FA = ifft2(TA);
    u = Vr - abs(FA).^2;
    Psi_x = exp(-p2*dt*u).*FA;
    
    % evolving in spatial space
    Psi_k = fft2(Psi_x);
    TA = A.*Psi_k;
    
    Psi = ifft2(TA);
    if mod(ind,1000) == 0 
        Na = sum(sum(abs(Psi).^2))*dx*dy;
        r2_mean = sum(sum(r.^2.*abs(Psi).^2))*dx*dy/Na;
        X_mean = sum(sum(X.*abs(Psi).^2))*dx*dy/Na;
        Y_mean = sum(sum(Y.*abs(Psi).^2))*dx*dy/Na;
        delta_r(indm) = r2_mean - (X_mean^2+Y_mean^2);
        indm = indm + 1;
    end
end