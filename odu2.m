function f = odu2(z, P, Lambda, sigma, psi, N, n_sum, ph_const, Pin, r_edf, N0)
%% итоговые скоростные уравнения (общий случай)
% Переменные
% r   - массив координат радиуса волокна
% N_r - размер массива координат радиуса волокна 
% n1  - населенность нижнего уровня
% n2  - населенность верхнего уровня
% f   - значение производных в скоростных уравнениях

% Функции
% population(P,...,r_edf) - расчет населенностей энергетических уровней в зависимости от мощности и длины волны

r  = 0 : r_edf / N0.r : r_edf;                                         % массив координат радиуса волокна
N_r = size(r,2);                                                       % размер массива координат радиуса волокна

[n1, n2]   = population(P, Lambda, sigma, psi, N, ph_const, Pin);      % расчет населенностей уровней

%уравнения для сигнала (dP_s/dz)
f(1:N.s,1) = 2 * pi * n_sum * trapz(r,(repmat(sigma.es',1,N_r) .* repmat(n2,N.s,1)...
    - repmat(sigma.as',1,N_r) .* repmat(n1,N.s,1)).* repmat(P(1:N.s,1),1,N_r) .* psi.ns .* repmat(r,N.s,1),2);

% уравнения для попутного ASE (dP_ase/dz)
f(N.s+1:N.s+N.ase,1)  = 2 * pi * n_sum * (trapz(r,(repmat(sigma.ease',1,N_r) .* repmat(n2,N.ase,1)...
    .*(repmat(P(N.s+1:N.s+N.ase,1),1,N_r) + repmat(2 * ph_const.h * ph_const.c^2 * 0.1 * 10^(-9) ./...
    Lambda.ase'.^3,1,N_r)) - repmat(sigma.aase',1,N_r) .* repmat(n1,N.ase,1) .* repmat(P(N.s+1:N.s+N.ase,1),1,N_r))...
    .* psi.nase .* repmat(r,N.ase,1),2));

% уравнения для попутной накачки (dP_p/dz)
f(N.s+N.ase+1:N.s+N.ase+N.pf,1)  = 2 * pi * n_sum * (trapz(r,(repmat(sigma.epf',1,N_r) .* repmat(n2,N.pf,1)...
    - repmat(sigma.apf',1,N_r) .* repmat(n1,N.pf,1)) .* repmat(P(N.s+N.ase+1:N.s+N.ase+N.pf,1),1,N_r) .*...
    psi.npf .* repmat(r,N.pf,1),2));

% допольнительные уравнения при встречной накачке
if isempty(Pin.pb) == 0
    % уравнения для встречного ASE (dP_ase/dz)
    f(N.s+N.ase+N.pf+1:N.s+N.ase+N.pf+N.ase,1) = -2 * pi * n_sum * (trapz(r,(repmat(sigma.ease',1,N_r) .*...
        repmat(n2,N.ase,1) .* (repmat(P(N.s+N.ase+N.pf+1:N.s+N.ase+N.pf+N.ase,1),1,N_r) + ...
        repmat(2 * ph_const.h * ph_const.c^2 * 0.1 * 10^(-9) ./ Lambda.ase'.^3,1,N_r)) - repmat(sigma.aase',1,N_r) .*...
        repmat(n1,N.ase,1) .* repmat(P(N.s+N.ase+N.pf+1:N.s+N.ase+N.pf+N.ase,1),1,N_r)) .* psi.nase .* repmat(r,N.ase,1),2));
    
    % уравнение для встречной накачки (dP_p/dz)
    f(N.s+N.ase+N.pf+N.ase+1:N.s+N.ase+N.pf+N.ase+N.pb,1) = -2 * pi * n_sum * (trapz(r,(repmat(sigma.epb',1,N_r) .*...
        repmat(n2,N.pb,1) - repmat(sigma.apb',1,N_r) .* repmat(n1,N.pb,1)) .*...
        repmat(P(N.s+N.ase+N.pf+N.ase+1:N.s+N.ase+N.pf+N.ase+N.pb,1),1,N_r) .* psi.npb .* repmat(r,N.pb,1),2));
end
end