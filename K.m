f = dir('th*');
mu = 0.1666;
G = [];
V = [];
for n=1:3
    cd(f(n).name);
    G = [G importdata('param.res')'];
    V = [V importdata('flux.res').data(end,2:end)'];
    dims  = importdata('dimensions.res');
    Gamma = h5read('final.h5','/Gamma');
    Por   = sum(Gamma)/length(Gamma);
    cd ..;
end

k = mu*V*inv(G);
k12 = 0.5*(k(1,2)+k(2,1));
k13 = 0.5*(k(1,3)+k(3,1));
k23 = 0.5*(k(2,3)+k(3,2));
k(2,1) = k12;
k(1,2) = k12;
k(3,1) = k13;
k(1,3) = k13;
k(3,2) = k23;
k(2,3) = k23;
[V,D] = eig(k);
kv1   = real(D(1,1))*real(V(:,1));
kv2   = real(D(2,2))*real(V(:,2));
kv3   = real(D(3,3))*real(V(:,3));
kv4   = -kv1;
kv5   = -kv2;
kv6   = -kv3;

