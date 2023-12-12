Lx = 8;
Ly = 8;

dx = 0.05;
dy = 0.05;
x = -Lx:dx:Lx;
y = -Ly:dy:Ly;

[X,Y] = meshgrid(x,y);
[phi,r] = cart2pol(X,Y);

Vr = 0.5*r.^2;

n = 12;
l = -4;
Psi = genEigenMode(n,l,r,phi);

[Fx,Fy] = gradient(Psi,dx,dy);
[Fxx,~] = gradient(Fx,dx,dy);
[~,Fyy] = gradient(Fy,dx,dy);

% compute the angular momentum, should be l, if correct
px = -1i*Fx;
py = -1i*Fy;
Lz = X.*py - Y.*px;
disp(['The numerical angular momentum is ',num2str(real(sum(sum(conj(Psi).*Lz))*dx*dy))])
disp(['The exact angular momentum is ',num2str(l)])

% compute the eigen-energy, should be n+1, if correct
H = -0.5*Fxx - 0.5*Fyy + Vr.*Psi;
E = sum(sum(conj(Psi).*H))*dx*dy;
disp(['The numerical eigen-energy is ', num2str(real(E))]);
disp(['The exact eigen-energy is ',num2str(n+1)])