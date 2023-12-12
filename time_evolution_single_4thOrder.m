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

% Suzuki-Trotter coefficient
s = 1/(2-2^(1/3));
p1 = s/2;
p2 = s;
p3 = (1-s)/2;
p4 = (1-2*s);
p5 = (1-s)/2;
p6 = s;
p7 = s/2;


load basis.mat;
M = size(basis);
M = M(1);

NN = 1000000;
MM = 1000;
w0 = zeros(1,NN/MM);
w1 = zeros(1,NN/MM);
w2 = zeros(1,NN/MM);
w3 = zeros(1,NN/MM);
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

N_mode = 12;
Z = 0:N_mode;
alpha = 2;

c = exp(-abs(alpha)^2/2).*alpha.^Z./sqrt(factorial(Z));

Psi0 = 0;

for ind = 1:N_mode+1
    Psi0 = Psi0 + c(ind)*squeeze(basis(ind,:,:));
end

% real time evolution, time step
dt = 1i*0.01;
% initial wave-function
Psi = Psi0;
Pn = zeros(1,NN/MM);
A = exp(-p7*T*dt);
B = exp(-p3*T*dt);
indm = 1;
for ind = 0:NN
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
    
    if mod(ind,MM) == 0
        ind
        w0(indm) = abs(sum(sum(conj(phi0).*Psi))*dx*dy)^2;
        w1(indm) = abs(sum(sum(conj(phi1).*Psi))*dx*dy)^2;
        w2(indm) = abs(sum(sum(conj(phi2).*Psi))*dx*dy)^2;
        w3(indm) = abs(sum(sum(conj(phi3).*Psi))*dx*dy)^2;
        Pn(indm) = sum(sum(abs(Psi).^2))*dx*dy;
        indm = indm + 1;        
    end
    
end
