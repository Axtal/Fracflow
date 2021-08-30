
KV1 = [];
KV2 = [];
KV3 = [];
KV4 = [];
KV5 = [];
KV6 = [];
X   = [];
Y   = [];
Z   = [];
P   = [];
idx = 1;
for nx0 = 100:200:300
    foldername = sprintf('nx0_%d',nx0);
    cd(foldername);
    K
    KV1 = [KV1 kv1];
    KV2 = [KV2 kv2];
    KV3 = [KV3 kv3];
    KV4 = [KV4 kv4];
    KV5 = [KV5 kv5];
    KV6 = [KV6 kv6];
    X   = [X   nx0];
    Y   = [Y   100];
    Z   = [Z   100];
    P   = [P   Por];
    cd ..
end

figure(1)
hold on;
axis square;
quiver3(X,Y,Z,KV1(1,:),KV1(2,:),KV1(3,:));
quiver3(X,Y,Z,KV2(1,:),KV2(2,:),KV2(3,:));
quiver3(X,Y,Z,KV3(1,:),KV3(2,:),KV3(3,:));
quiver3(X,Y,Z,KV4(1,:),KV4(2,:),KV4(3,:));
quiver3(X,Y,Z,KV5(1,:),KV5(2,:),KV5(3,:));
quiver3(X,Y,Z,KV6(1,:),KV6(2,:),KV6(3,:));
