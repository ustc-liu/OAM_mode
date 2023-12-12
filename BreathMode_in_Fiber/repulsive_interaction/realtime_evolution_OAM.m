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

kappa = 2;
Vr = 0.5*r.^2;
Vt = 0.5*kappa^2*r.^2;
% interaction strength

g = 1;

% Suzuki-Trotter coefficients
s = 1/(2-2^(1/3));
p1 = s/2;
p2 = s;
p3 = (1-s)/2;
p4 = (1-2*s);
p5 = (1-s)/2;
p6 = s;
p7 = s/2;

% imagnary-time evolution, time step
dt = 1i*0.001;
Nt1 = 200;
t1 = Nt1*(-1i*dt);
if kappa == 0
    t0 = pi/2 - atan(t1/2);
else
    t0 = pi/2 - atan((1+kappa^2)*tan(kappa*t1)/(2*kappa));
end
Nt2 = floor(t0/dt*1i);


% initial wave-function
load basis.mat;
N_mode = 1;
c = exp(1i*2*pi*rand(N_mode,1))/sqrt(N_mode);
Psi0 = 0;
Psi0 = Psi0 + c*squeeze(basis(2,:,:));
Psi = Psi0;

A = exp(-p7*T*dt);
B = exp(-p3*T*dt);
% number of T 
NT = 9;

delta_r = zeros(NT*(Nt1+Nt2),1);

% overlap of the wavefunctions;
Ft = zeros(NT*(Nt1+Nt2),1);

indm = 1;
for inda = 1:NT
    %% the first evolution duration
    for ind = 0:Nt1-1
        % evolving in momentum space
        Psi_k = fft2(Psi);
        TA = A.*Psi_k;
        
        % evolving in spatial space
        FA = ifft2(TA);
        % non-linear and potential part
        u = Vt + g*abs(FA).^2;
        Psi_x = exp(-p6*dt*u).*FA;
        
        % evolving in momentum space
        Psi_k = fft2(Psi_x);
        TA = B.*Psi_k;
        
        % evolving in momentum space
        FA = ifft2(TA);
        u = Vt + g*abs(FA).^2;
        Psi_x = exp(-p4*dt*u).*FA;
        
        % evolving in spatial space
        Psi_k = fft2(Psi_x);
        TA = B.*Psi_k;
        
        % evolving in momentum space
        FA = ifft2(TA);
        u = Vt + g*abs(FA).^2;
        Psi_x = exp(-p2*dt*u).*FA;
        
        % evolving in spatial space
        Psi_k = fft2(Psi_x);
        TA = A.*Psi_k;
        
        Psi = ifft2(TA);
        
        Na = sum(sum(abs(Psi).^2))*dx*dy;
        r2_mean = sum(sum(r.^2.*abs(Psi).^2))*dx*dy/Na;
        X_mean = sum(sum(X.*abs(Psi).^2))*dx*dy/Na;
        Y_mean = sum(sum(Y.*abs(Psi).^2))*dx*dy/Na;
        delta_r(indm) = sqrt(r2_mean - (X_mean^2+Y_mean^2));
        Ft(indm) = abs(sum(sum(conj(Psi0).*Psi))*dx*dy)^2;
        indm = indm + 1;
        
    end
    %% the second evolution duration
    for ind = 0:Nt2-1
        % evolving in momentum space
        Psi_k = fft2(Psi);
        TA = A.*Psi_k;
        
        % evolving in spatial space
        FA = ifft2(TA);
        % non-linear and potential part
        u = Vr + g*abs(FA).^2;
        Psi_x = exp(-p6*dt*u).*FA;
        
        % evolving in momentum space
        Psi_k = fft2(Psi_x);
        TA = B.*Psi_k;
        
        % evolving in momentum space
        FA = ifft2(TA);
        u = Vr + g*abs(FA).^2;
        Psi_x = exp(-p4*dt*u).*FA;
        
        % evolving in spatial space
        Psi_k = fft2(Psi_x);
        TA = B.*Psi_k;
        
        % evolving in momentum space
        FA = ifft2(TA);
        u = Vr + g*abs(FA).^2;
        Psi_x = exp(-p2*dt*u).*FA;
        
        % evolving in spatial space
        Psi_k = fft2(Psi_x);
        TA = A.*Psi_k;
        
        Psi = ifft2(TA);
        
        Na = sum(sum(abs(Psi).^2))*dx*dy;
        r2_mean = sum(sum(r.^2.*abs(Psi).^2))*dx*dy/Na;
        X_mean = sum(sum(X.*abs(Psi).^2))*dx*dy/Na;
        Y_mean = sum(sum(Y.*abs(Psi).^2))*dx*dy/Na;
        delta_r(indm) = sqrt(r2_mean - (X_mean^2+Y_mean^2));
        Ft(indm) = abs(sum(sum(conj(Psi0).*Psi))*dx*dy)^2;
        indm = indm + 1;
    end
end
subplot(2,1,1)
plot((-1i*dt)*(0:NT*(Nt1+Nt2)-1),Ft)
hold on
plot((-1i*dt)*2*(Nt1+Nt2)*ones(100,1),linspace(min(Ft),1,100),'--r')
plot((-1i*dt)*4*(Nt1+Nt2)*ones(100,1),linspace(min(Ft),1,100),'--r')

subplot(2,1,2)
plot((-1i*dt)*(0:NT*(Nt1+Nt2)-1),delta_r)
hold on
plot((-1i*dt)*2*(Nt1+Nt2)*ones(100,1),linspace(min(delta_r),max(delta_r),100),'--r')
plot((-1i*dt)*4*(Nt1+Nt2)*ones(100,1),linspace(min(delta_r),max(delta_r),100),'--r')