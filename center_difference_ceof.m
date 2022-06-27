clc
clear

dx0 = 0;
order = 8;
np = order + 1; % Number of source points
nt = np; % Number of polynomial terms
ic = order / 2 + 1;

x = zeros(np,1);
for ip = 1:np
    x(ip) = ip - ic;
end

poly = zeros(np,1);
for ip = 1:np
    for it = 1:nt
        poly(ip,it) = x(ip).^(it-1);
    end
end

dpoly1 = zeros(1,nt);
for it = 1:nt
    dpoly1(1,it) = max(it-1,0) * dx0.^max([it-2,0]);
end

coef = dpoly1 / poly;