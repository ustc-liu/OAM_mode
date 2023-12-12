xmax = 8;
ymax = 8;
Nx = 2^7;
Ny = 2^7;
dx = 2*xmax/Nx;
dy = 2*ymax/Ny;
x = -xmax:dx:xmax-dx;
y = -ymax:dy:ymax-dy;

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

load basis.mat;
M = size(basis);
M = M(1);

NN = 100000;

dt = 1i*0.001;

N_norm = 30;
%c = sqrt(N_norm)*exp(1i*rand(M,1))./sqrt(M);
c = sqrt(N_norm);
Psi0 = 0;
for ind = 1:1
    Psi0 = Psi0 + c(ind)*squeeze(basis(ind,:,:));
end

% Psi0 = exp(-r.^2);
% Psi0 = Psi0/sqrt(sum(sum(abs(Psi0).^2))*dx*dy)*sqrt(N_norm);

%% evolution
% initial wave-function
Psi = Psi0;
A = exp(-0.5*T*dt);
delta_r = zeros(1,NN/1000);
indz = 1;
for ind = 0:NN-1
    
    if mod(ind,1000) == 0
        Nt = sum(sum(abs(Psi).^2))*dx*dy;
        r2_mean = sum(sum(r.^2.*abs(Psi).^2))*dx*dy/Nt;
        x_mean = sum(sum(X.*abs(Psi).^2))*dx*dy/Nt;
        y_mean = sum(sum(Y.*abs(Psi).^2))*dx*dy/Nt;
        r_mean2 = x_mean^2 + y_mean^2;
        delta_r(indz) = r2_mean - r_mean2;
        indz = indz + 1;
    end
    % evolving in momentum space
    TA = A.*fft2(Psi);
    
    % evolving in spatial space
    FA = ifft2(TA);
    % non-linear and potential part
    u = Vr - abs(FA).^2;
    Psi_x = exp(-dt*u).*FA;
    
    % evolving in momentum space
    TA = A.*fft2(Psi_x);
    
    Psi = ifft2(TA);
    
end

% Psi0 = sqrt(N_norm)*Psi0;
% Psi = sqrt(N_norm)*Psi;

Psik0 = fft2(Psi0)*dx*dy/(2*pi);
%Nk0 = sum(sum(abs(Psik0).^2))*dkx*dky;
%Psik0 = Psik0./sqrt(Nk0)*sqrt(N_norm);
%
%
% Nk = sum(sum(abs(TA).^2))*dkx*dky;
% TA = TA./sqrt(Nk)*sqrt(N_norm);
Psik = fft2(Psi)*dx*dy/(2*pi);
%Nk = sum(sum(abs(Psik).^2))*dkx*dky;
%Psik = Psik/sqrt(Nk)*sqrt(N_norm);

%% initial state energy
% initial kinetic energy
Ek = sum(sum(abs(Psik0).^2.*T))*dkx*dky;
% initial potential energy
Ep = sum(sum(Vr.*abs(Psi0).^2))*dx*dy;
% initial interaction energy
Eint = -0.5*sum(sum(abs(Psi0).^4))*dx*dy;
E_start = Ek + Ep + Eint;

%% final state energy
Ek_f = sum(sum(abs(Psik).^2.*T))*dkx*dky;
% initial potential energy
Ep_f = sum(sum(Vr.*abs(Psi).^2))*dx*dy;
% initial interaction energy
Eint_f = -0.5*sum(sum(abs(Psi).^4))*dx*dy;
E_final = Ek_f + Ep_f + Eint_f;