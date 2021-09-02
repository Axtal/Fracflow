sol = h5read('Coral-Trial.h5','/t0/channel0');
dims = size(sol);

nx0 = 201;
ny0 = 201;
nz0 = 1;
nx  = 101;
ny  = 101;
nz  = 11;

ix = dims(3)-(nx0:nx0+nx-1);
iy = dims(2)-(ny0:ny0+ny-1);
iz = dims(1)-(nz0:nz0+nz-1);

solpl = sol(iz,iy,ix);

[x,y,z] = meshgrid(ix,iy,iz);

idx = find(solpl==1);

plot3(x(idx),y(idx),z(idx),'o');

