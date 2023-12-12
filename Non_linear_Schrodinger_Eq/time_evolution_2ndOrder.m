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

NN = 10000;

phi0 = squeeze(basis(1,:,:));
phi1 = squeeze(basis(2,:,:));
phi2 = squeeze(basis(3,:,:));
phi3 = squeeze(basis(4,:,:));

% t0 = 10;
epsilon = [30,28,26,24,22,20,18,16,14,12,10,8];
E_indx = {1,[2,3],[4,5,6],[7,8,9,10],[11,12,13,14,15],[16,17,18,19,20,21],[22,23,24,25,26,27,28],[29,30,31,32,33,34,35,36]...
    [37,38,39,40,41,42,43,44,45],[46,47,48,49,50,51,52,53,54,55],[56,57,58,59,60,61,62,63,64,65,66],[67,68,69,70,71,72,73,74,...
    75,76,77,78]};

% new_basis = zeros(M,Nx,Ny);
%
% for ind = 1:12
%     en = epsilon(ind);
%     theta = exp(1i*en*t0);
%     for indx = E_indx{ind}
%         new_basis(indx,:,:) = basis(indx,:,:)*theta;
%     end
% end
% real time evolution, time step
dt = 1i*0.0001;
w0 = zeros(10,NN);
w2 = zeros(10,NN);
Pn = zeros(10,NN);

parfor indz = 1:10
    c = exp(1i*rand(1,M)*2*pi)/sqrt(M);
    c(1) = sqrt(1/6)*exp(1i*rand*2*pi);
    c(2) = sqrt(1/7)*exp(1i*rand*2*pi);
    c(3) = sqrt(1/8)*exp(1i*rand*2*pi);
    c(4) = sqrt(1/9)*exp(1i*rand*2*pi);
    c(5) = sqrt(1/10)*exp(1i*rand*2*pi);
    c(6) = sqrt(1/10)*exp(1i*rand*2*pi);
    c = c/sqrt(sum(abs(c).^2));
    
    Psi0 = 0;
    
    for ind = 1:M
        Psi0 = Psi0 + c(ind)*squeeze(basis(ind,:,:));
    end
    
    
    % initial wave-function
    Psi = Psi0;
    A = exp(-0.5*T*dt);
    for ind = 1:NN
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
        
        w0(indz,ind) = abs(sum(sum(conj(phi0).*Psi))*dx*dy)^2;
        w2(indz,ind) = abs(sum(sum(conj(phi2).*Psi))*dx*dy)^2;
        Pn(indz,ind) = sum(sum(abs(Psi).^2))*dx*dy;
    end
end