function f = odu2(z, P, wl, sigma, psi, N, n_sum, ph_const, P_in, r_edf)
%% итоговые скоростные уравнения (общий случай)
% Переменные
% r   - массив координат радиуса волокна
% N_r - размер массива координат радиуса волокна 
% n1  - населенность нижнего уровня
% n2  - населенность верхнего уровня
% f   - значение производных в скоростных уравнениях

% Функции
% population(P,...,ph_const) - расчет населенностей энергетических уровней в зависимости от мощности и длины волны

r = 0 : 10^(-7) : 10^(-5);                                                 % массив координат радиуса волокна
N_r = size(r,2);                                                           % размер массива координат радиуса волокна

[n1, n2]   = population(P, wl, sigma, psi, N, ph_const, P_in, r_edf);                   % расчет населенностей уровней


%уравнения для сигнала (dP_s/dz)
f(1:N.S,1) = 2 * pi * n_sum * trapz(r,(repmat(sigma.ES',1,N_r) .* repmat(n2,N.S,1)...
    - repmat(sigma.AS',1,N_r) .* repmat(n1,N.S,1)).* repmat(P(1: N.S, 1),1,N_r) .* psi.NS .* repmat(r,N.S,1),2);

% уравнения для попутного ASE (dP_ase/dz)
f(N.S+1:N.S+N.ASE,1)  = 2 * pi * n_sum * (trapz(r,(repmat(sigma.EASE',1,N_r) .* repmat(n2,N.ASE,1)...
    .*(repmat(P(N.S+1: N.S+N.ASE, 1),1,N_r) + repmat(2 * ph_const.h * ph_const.c^2 * 0.1 * 10^18 ./...
    wl.ASE'.^3,1, N_r)) - repmat(sigma.AASE',1,N_r) .* repmat(n1,N.ASE,1) .* repmat(P(N.S+1: N.S+N.ASE, 1),1,N_r))...
    .* psi.NASE .* repmat(r,N.ASE,1),2));

% уравнения для попутной накачки (dP_p/dz)
f(N.S+N.ASE+1:N.S+N.ASE+N.PF,1)  = 2 * pi * n_sum * (trapz(r,(repmat(sigma.EPF',1,N_r) .* repmat(n2,N.PF,1)...
    - repmat(sigma.APF',1,N_r) .* repmat(n1,N.PF,1)) .* repmat(P(N.S+N.ASE+1: N.S+N.ASE+N.PF, 1),1,N_r) .*...
    psi.NPF .* repmat(r,N.PF,1),2));

% допольнительные уравнения при встречной накачке
if isempty(P_in.PB) == 0
    % уравнения для встречного ASE (dP_ase/dz)
    f(N.S+N.ASE+N.PF+1:N.S+N.ASE+N.PF+N.ASE,1) = -2 * pi * n_sum * (trapz(r,(repmat(sigma.EASE',1,N_r) .*...
        repmat(n2,N.ASE,1) .* (repmat(P(N.S+N.ASE+N.PF+1: N.S+N.ASE+N.PF+N.ASE, 1),1,N_r) + ...
        repmat(2 * ph_const.h * ph_const.c^2 * 0.1 * 10^18 ./ wl.ASE'.^3,1,N_r)) - repmat(sigma.AASE',1,N_r) .*...
        repmat(n1,N.ASE,1) .* repmat(P(N.S+N.ASE+N.PF+1: N.S+N.ASE+N.PF+N.ASE, 1),1,N_r)) .* psi.NASE .* repmat(r,N.ASE,1),2));
    
    % уравнение для встречной накачки (dP_p/dz)
    f(N.S+N.ASE+N.PF+N.ASE+1:N.S+N.ASE+N.PF+N.ASE+N.PB,1) = -2 * pi * n_sum * (trapz(r,(repmat(sigma.EPB',1,N_r) .*...
        repmat(n2,N.PB,1) - repmat(sigma.APB',1,N_r) .* repmat(n1,N.PB,1)) .*...
        repmat(P(N.S+N.ASE+N.PF+N.ASE+1: N.S+N.ASE+N.PF+N.ASE+N.PB, 1),1,N_r) .* psi.NPB .* repmat(r,N.PB,1),2));
end
end