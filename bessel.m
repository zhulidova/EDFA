%% –асчет волновых функций распределени€ интенсивности в волокне
% зависимости от радиуса волокна дл€ каждой длины волны
function [psi, w_edf] = bessel(r_edf, wl, N, NA_edf)
% ѕеременные
% V_edf         - структура нормированных частот дл€ сигнала (V_edf.S), накачки (V_edf.PF и
%V_edf.PB) и ASE (V_edf.ASE)
% U_edf и W_edf - структура параметров дл€ фукций Ѕессел€ в зависимости от длины волны
% w_edf         - структура модовых площадей дл€ сигнала (w_edf.S), накачки (w_edf.PF и w_edf.PB) и ASE (w_edf.ASE)

% расчет нормированных частот
V_edf.S                = 2 * pi * r_edf * NA_edf * 10^9 ./ wl.S;    % нормированна€ частота сигнала
V_edf.PF               = 2 * pi * r_edf * NA_edf * 10^9 ./ wl.PF;   % нормированна€ частота попутной накачки
V_edf.PB               = 2 * pi * r_edf * NA_edf * 10^9 ./ wl.PB;   % нормированна€ частота встречной накачки
V_edf.ASE              = 2 * pi * r_edf * NA_edf * 10^9 ./ wl.ASE;  % нормированна€ частота ASE

% параметры дл€ расчета функций Ѕессел€
U_edf.S                = (1 + sqrt(2)) .* V_edf.S ./ (1 + (4 + V_edf.S.^4).^0.25);
U_edf.PF               = (1 + sqrt(2)) .* V_edf.PF ./ (1 + (4 + V_edf.PF.^4).^0.25);
U_edf.PB               = (1 + sqrt(2)) .* V_edf.PB ./ (1 + (4 + V_edf.PB.^4).^0.25);
U_edf.ASE              = (1 + sqrt(2)) .* V_edf.ASE ./ (1 + (4 + V_edf.ASE.^4).^0.25);
W_edf.S                = (V_edf.S.^2 - U_edf.S.^2).^0.5;
W_edf.PF               = (V_edf.PF.^2 - U_edf.PF.^2).^0.5;
W_edf.PB               = (V_edf.PB.^2 - U_edf.PB.^2).^0.5;
W_edf.ASE              = (V_edf.ASE.^2 - U_edf.ASE.^2).^0.5;

% модовые площади
w_edf.S                = r_edf .* V_edf.S .* besselk(1, W_edf.S)...    % модова€ площадь (сигнал)
    .* besselj(0, U_edf.S) ./ (U_edf.S .* besselk(0, W_edf.S)); 
w_edf.PF               = r_edf .* V_edf.PF .* besselk(1, W_edf.PF)...  % модова€ площадь (попутна€ накачка)
    .* besselj(0, U_edf.PF) ./ (U_edf.PF .* besselk(0, W_edf.PF));
w_edf.PB               = r_edf .* V_edf.PB .* besselk(1, W_edf.PB)...  % модова€ площадь (встречна€ накачка)
    .* besselj(0, U_edf.PB) ./ (U_edf.PB .* besselk(0, W_edf.PB));
w_edf.ASE              = r_edf .* V_edf.ASE .* besselk(1, W_edf.ASE)...% модова€ площадь (ASE)
    .* besselj(0, U_edf.ASE) ./ (U_edf.ASE .* besselk(0, W_edf.ASE));

% расчет зависимости функций распределени€ интенсивности от рассто€ни€ от центра волокна
% дл€  r < r_edf
r_less_a = 0 : 10^(-7) : r_edf;
n_r_la   = size(r_less_a,2);                                           % размер массива

psi.S(1: N.S,1: n_r_la)     = (besselj(0, (repmat(U_edf.S',1,n_r_la)...  
    .* repmat(r_less_a, N.S,1)) / r_edf)).^2; 
psi.PF(1: N.PF,1: n_r_la)   = (besselj(0, (repmat(U_edf.PF',1,n_r_la)...
    .* repmat(r_less_a, N.PF,1)) / r_edf)).^2;
psi.PB(1: N.PB,1: n_r_la)   = (besselj(0, (repmat(U_edf.PB',1,n_r_la)...
    .* repmat(r_less_a, N.PB,1)) / r_edf)).^2;
psi.ASE(1: N.ASE,1: n_r_la) = (besselj(0, (repmat(U_edf.ASE',1,n_r_la)...
    .* repmat(r_less_a, N.ASE,1)) / r_edf)).^2;
psi.NS             = psi.S ./ pi ./ w_edf.S'.^2;                       % нормированна€ функци€ распределени€ (сигнал)
psi.NPF            = psi.PF ./ pi ./ w_edf.PF'.^2;                     % нормированна€ функци€ распределени€ (попутна€ накачка)
psi.NPB            = psi.PB ./ pi ./ w_edf.PB'.^2;                     % нормированна€ функци€ распределени€ (встречна€ накачка)
psi.NASE           = psi.ASE ./ pi ./ w_edf.ASE'.^2;                   % нормированна€ функци€ распределени€ (ASE)

% дл€  r>a
r_more_a = r_edf: 10^(-7): 10^(-5);
n_r_ma   = size(r_more_a,2);                                          % размер массива

psi.S(1: N.S,n_r_la+1: n_r_la+n_r_ma)     = (besselk(0,(repmat(W_edf.S',1,n_r_ma) .*...
    repmat(r_more_a, N.S,1)) / r_edf)).^2 .* besselj(0, repmat(U_edf.S',1,n_r_ma)).^2 ./...
    besselk(0, repmat(W_edf.S',1,n_r_ma)).^2;  
psi.PF(1: N.PF,n_r_la+1: n_r_la+n_r_ma)   = (besselk(0,(repmat(W_edf.PF',1,n_r_ma)...
    .* repmat(r_more_a, N.PF,1)) / r_edf)).^2 .* besselj(0, repmat(U_edf.PF',1,n_r_ma)).^2 ./...
    besselk(0, repmat(W_edf.PF',1,n_r_ma)).^2;
psi.PB(1: N.PB,n_r_la+1: n_r_la+n_r_ma)   = (besselk(0,(repmat(W_edf.PB',1,n_r_ma)...
    .* repmat(r_more_a, N.PB,1)) / r_edf)).^2 .* besselj(0, repmat(U_edf.PB',1,n_r_ma)).^2 ./...
    besselk(0, repmat(W_edf.PB',1,n_r_ma)).^2;
psi.ASE(1: N.ASE,n_r_la+1: n_r_la+n_r_ma) = (besselk(0,(repmat(W_edf.ASE',1,n_r_ma)...
    .* repmat(r_more_a, N.ASE,1)) / r_edf)).^2 .* besselj(0, repmat(U_edf.ASE',1,n_r_ma)).^2 ./...
    besselk(0, repmat(W_edf.ASE',1,n_r_ma)).^2;
psi.NS              = psi.S ./ pi ./ w_edf.S'.^2;    % нормированна€ функци€ распределени€ (сигнал)
psi.NPF             = psi.PF ./ pi ./ w_edf.PF'.^2;  % нормированна€ функци€ распределени€ (попутна€ накачка)
psi.NPB             = psi.PB ./ pi ./ w_edf.PB'.^2;  % нормированна€ функци€ распределени€ (встречна€ накачка)
psi.NASE            = psi.ASE ./ pi ./ w_edf.ASE'.^2;% нормированна€ функци€ распределени€ (ASE)

end