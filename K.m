f = dir('th*');
mu = 0.1666;
G = [];
V = [];
for n=1:length(f)
    cd(f(n).name);
    G = [G importdata('param.res')'];
    V = [V importdata('flux.res').data(100,2:end)'];
    cd ..;
end

k = mu*V*inv(G);