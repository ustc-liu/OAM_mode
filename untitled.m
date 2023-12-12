load basis.mat;
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
[M,~,~] = size(basis);
dV = zeros(1000,Nx,Ny);
for ind = 1:1000
    dn = randn(M,1);
    for indx = 1:M
        dV(ind,:,:) = dV(ind,:,:) + dn(indx)*abs(basis(indx,:,:)).^2;
    end
end