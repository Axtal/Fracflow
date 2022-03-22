sol = h5read('Coral-Trial.h5','/t0/channel0');
dims = size(sol);



nx0 = 0;
ny0 = 0;
nz0 = 100;
nx  = 600;
ny  = 600;
nz  = 1;
step = 1;

ixpt = dims(1)-(nx0:nx0+nx-1);
iypt = dims(2)-(ny0:ny0+ny-1);
izpt = dims(3)-(nz0:nz0+nz-1);
ix = (nx0:nx0+nx-1);
iy = (ny0:ny0+ny-1);
iz = (nz0:nz0+nz-1);

solpt = rot90(flip(sol(ixpt,iypt,izpt),1),-1);

[x,y,z] = meshgrid(ix,iy,iz);

idx=find(solpt==1);

xp = x(idx(1:step:end));
yp = y(idx(1:step:end));
zp = z(idx(1:step:end));

figure(1)
plot3(xp,yp,zp,'.'); %use this for z planes
hold on;
axis square;
xlabel('x');
ylabel('y');
zlabel('z');



