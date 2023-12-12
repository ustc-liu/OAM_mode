xmax = 8;
ymax = 8;
Nx = 2^8;
Ny = 2^8;
dx = 2*xmax/Nx;
dy = 2*ymax/Ny;
x = -xmax:dx:xmax-dx;
y = -ymax:dy:ymax-dy;

% space-spatial meshgrid
[X,Y] = meshgrid(x,y);
[phi,r] = cart2pol(X,Y);

L = {0,[-1,1],[-2,0,2],[-3,-1,1,3],[-4,-2,0,2,4],[-5,-3,-1,1,3,5],[-6,-4,-2,0,2,4,6],[-7,-5,-3,-1,1,3,5,7],[-8,-6,-4,-2,0,2,4,6,8]};
M = 0;
for ind = 1:length(L)
    M = M + length(L{ind});
end
basis = zeros(M,Nx,Ny);
ind = 1;
for indx = 1:length(L)
    for indy = 1:length(L{indx})
        basis(ind,:,:) = genEigenMode(indx-1,L{indx}(indy),r,phi);
        ind = ind + 1;
    end
end
save('basis.mat','basis');
