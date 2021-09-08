sol = h5read('Coral-Trial.h5','/t0/channel0');
dims = size(sol);

nx0 = 200;
ny0 = 300;
nz0 = 50;
nx  = 300;
ny  = 300;
nz  = 1;
step = 1;

ixpt = dims(1)-(nx0:nx0+nx-1);
iypt = dims(2)-(ny0:ny0+ny-1);
izpt = dims(3)-(nz0:nz0+nz-1);
ix = (nx0:nx0+nx-1);
iy = (ny0:ny0+ny-1);
iz = (nz0:nz0+nz-1);

solpt = sol(ixpt,iypt,izpt);

[x,y,z] = meshgrid(ix,iy,iz);

idx=find(solpt==1);

% xp = -y(idx)+min(y(idx))+max(y(idx));
% yp = -x(idx)+min(x(idx))+max(x(idx));
xt = x(idx(1:step:end));
yt = y(idx(1:step:end));
zp = z(idx(1:step:end));

xp = nx0 +(yt-ny0);
yp = ny0 +(xt-nx0);

figure(1)
plot3(xp,yp,zp,'o'); %use this for z planes
%plot3(xt,yt,zp,'o'); %use this for y planes
hold on;
axis square;
xlabel('x');
ylabel('y');
zlabel('z');



