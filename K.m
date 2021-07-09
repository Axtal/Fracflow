f = dir('th*');
mu = 0.1666;
G = [];
V = [];
for n=4:6
    cd(f(n).name);
    G = [G importdata('param.res')'];
    V = [V importdata('flux.res').data(end,2:end)'];
    cd ..;
end

k = mu*V*inv(G);