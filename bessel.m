%% Расчет волновых функций распределения интенсивности в волокне
% зависимости от радиуса волокна для каждой длины волны
function [psi, w_edf] = bessel(r_edf, Lambda, N, NA_edf, N0)
% Переменные
% V_edf - структура нормированных частот для сигнала (V_edf.s), накачки (V_edf.pf и V_edf.pb) и ASE (V_edf.ase)
% U_edf и W_edf - структура параметров для фукций Бесселя в зависимости от длины волны
% w_edf - структура модовых радиусов для сигнала (w_edf.s), накачки (w_edf.pf и w_edf.PB) и ASE (w_edf.ASE)

% расчет нормированных частот
V_edf.s                = 2 * pi * r_edf * NA_edf ./ Lambda.s;    % нормированная частота сигнала
V_edf.pf               = 2 * pi * r_edf * NA_edf ./ Lambda.pf;   % нормированная частота попутной накачки
V_edf.ase              = 2 * pi * r_edf * NA_edf ./ Lambda.ase;  % нормированная частота ASE

% параметры для расчета функций Бесселя
U_edf.s                = (1 + sqrt(2)) .* V_edf.s ./ (1 + (4 + V_edf.s.^4).^0.25);
U_edf.pf               = (1 + sqrt(2)) .* V_edf.pf ./ (1 + (4 + V_edf.pf.^4).^0.25);
U_edf.ase              = (1 + sqrt(2)) .* V_edf.ase ./ (1 + (4 + V_edf.ase.^4).^0.25);
W_edf.s                = (V_edf.s.^2 - U_edf.s.^2).^0.5;
W_edf.pf               = (V_edf.pf.^2 - U_edf.pf.^2).^0.5;
W_edf.ase              = (V_edf.ase.^2 - U_edf.ase.^2).^0.5;

% модовые радиусы
w_edf.s                = 1.2*r_edf .* V_edf.s .* besselk(1, W_edf.s)...    % модовый радиус (сигнал)
    .* besselj(0, U_edf.s) ./ (U_edf.s .* besselk(0, W_edf.s)); 
w_edf.pf               = 1.2*r_edf .* V_edf.pf .* besselk(1, W_edf.pf)...  % модовый радиус (попутная накачка)
    .* besselj(0, U_edf.pf) ./ (U_edf.pf .* besselk(0, W_edf.pf));
w_edf.ase              = 1.2*r_edf .* V_edf.ase .* besselk(1, W_edf.ase)...% модовый радиус (ASE)
    .* besselj(0, U_edf.ase) ./ (U_edf.ase .* besselk(0, W_edf.ase));

% расчет зависимости функций распределения интенсивности от расстояния от центра волокна
% для  r < r_edf
r  = 0 : r_edf / N0 : r_edf;
n_r   = size(r,2);                           % размер массива

psi.s(1: N.s,1: n_r)     = (besselj(0, (repmat(U_edf.s',1,n_r) .* repmat(r, N.s,1)) / r_edf)).^2; 
psi.pf(1: N.pf,1: n_r)   = (besselj(0, (repmat(U_edf.pf',1,n_r) .* repmat(r, N.pf,1)) / r_edf)).^2;
psi.ase(1: N.ase,1: n_r) = (besselj(0, (repmat(U_edf.ase',1,n_r) .* repmat(r, N.ase,1)) / r_edf)).^2;

psi.ns    = psi.s ./ pi ./ w_edf.s'.^2;      % нормированная функция распределения (сигнал)
psi.npf   = psi.pf ./ pi ./ w_edf.pf'.^2;    % нормированная функция распределения (попутная накачка)
psi.nase  = psi.ase ./ pi ./ w_edf.ase'.^2;  % нормированная функция распределения (ASE)

% расчет параметров для встречной накачки, если она есть
if isempty(Lambda.pb) == 0
    V_edf.pb = 2 * pi * r_edf * NA_edf ./ Lambda.pb;   % нормированная частота встречной накачки
    U_edf.pb = (1 + sqrt(2)) .* V_edf.pb ./ (1 + (4 + V_edf.pb.^4).^0.25);
    W_edf.pb = (V_edf.pb.^2 - U_edf.pb.^2).^0.5;
    w_edf.pb = r_edf .* V_edf.pb .* besselk(1, W_edf.pb)...  % модовый радиус (встречная накачка)
    .* besselj(0, U_edf.pb) ./ (U_edf.pb .* besselk(0, W_edf.pb));
    psi.pb(1: N.pb,1: n_r)   = (besselj(0, (repmat(U_edf.pb',1,n_r) .* repmat(r, N.pb,1)) / r_edf)).^2;
    psi.npb            = psi.pb ./ pi ./ w_edf.pb'.^2;       % нормированная функция распределения (встречная накачка)
end

end