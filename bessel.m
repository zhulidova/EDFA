%% –асчет волновых функций распределени€ интенсивности в волокне
% зависимости от радиуса волокна дл€ каждой длины волны
function [psi, w_edf] = bessel(r_edf, Lambda, N, NA_edf, N0)
% ѕеременные
% V_edf         - структура нормированных частот дл€ сигнала (V_edf.S), накачки (V_edf.PF и
%V_edf.PB) и ASE (V_edf.ASE)
% U_edf и W_edf - структура параметров дл€ фукций Ѕессел€ в зависимости от длины волны
% w_edf         - структура модовых площадей дл€ сигнала (w_edf.S), накачки (w_edf.PF и w_edf.PB) и ASE (w_edf.ASE)

% расчет нормированных частот
V_edf.s                = 2 * pi * r_edf * NA_edf ./ Lambda.s;    % нормированна€ частота сигнала
V_edf.pf               = 2 * pi * r_edf * NA_edf ./ Lambda.pf;   % нормированна€ частота попутной накачки
V_edf.ase              = 2 * pi * r_edf * NA_edf ./ Lambda.ase;  % нормированна€ частота ASE

% параметры дл€ расчета функций Ѕессел€
U_edf.s                = (1 + sqrt(2)) .* V_edf.s ./ (1 + (4 + V_edf.s.^4).^0.25);
U_edf.pf               = (1 + sqrt(2)) .* V_edf.pf ./ (1 + (4 + V_edf.pf.^4).^0.25);
U_edf.ase              = (1 + sqrt(2)) .* V_edf.ase ./ (1 + (4 + V_edf.ase.^4).^0.25);
W_edf.s                = (V_edf.s.^2 - U_edf.s.^2).^0.5;
W_edf.pf               = (V_edf.pf.^2 - U_edf.pf.^2).^0.5;
W_edf.ase              = (V_edf.ase.^2 - U_edf.ase.^2).^0.5;

% модовые площади
w_edf.s                = r_edf .* V_edf.s .* besselk(1, W_edf.s)...    % модова€ площадь (сигнал)
    .* besselj(0, U_edf.s) ./ (U_edf.s .* besselk(0, W_edf.s)); 
w_edf.pf               = r_edf .* V_edf.pf .* besselk(1, W_edf.pf)...  % модова€ площадь (попутна€ накачка)
    .* besselj(0, U_edf.pf) ./ (U_edf.pf .* besselk(0, W_edf.pf));
w_edf.ase              = r_edf .* V_edf.ase .* besselk(1, W_edf.ase)...% модова€ площадь (ASE)
    .* besselj(0, U_edf.ase) ./ (U_edf.ase .* besselk(0, W_edf.ase));

% расчет зависимости функций распределени€ интенсивности от рассто€ни€ от центра волокна
% дл€  r < r_edf
r  = 0 : r_edf / N0 : r_edf;
n_r   = size(r,2);                                           % размер массива

psi.s(1: N.s,1: n_r)     = (besselj(0, (repmat(U_edf.s',1,n_r)...  
    .* repmat(r, N.s,1)) / r_edf)).^2; 
psi.pf(1: N.pf,1: n_r)   = (besselj(0, (repmat(U_edf.pf',1,n_r)...
    .* repmat(r, N.pf,1)) / r_edf)).^2;
psi.ase(1: N.ase,1: n_r) = (besselj(0, (repmat(U_edf.ase',1,n_r)...
    .* repmat(r, N.ase,1)) / r_edf)).^2;
psi.ns             = psi.s ./ pi ./ w_edf.s'.^2;                       % нормированна€ функци€ распределени€ (сигнал)
psi.npf            = psi.pf ./ pi ./ w_edf.pf'.^2;                     % нормированна€ функци€ распределени€ (попутна€ накачка)
psi.nase           = psi.ase ./ pi ./ w_edf.ase'.^2;                   % нормированна€ функци€ распределени€ (ASE)

if isempty(Lambda.pb) == 0
    V_edf.pb               = 2 * pi * r_edf * NA_edf ./ Lambda.pb;   % нормированна€ частота встречной накачки
    U_edf.pb               = (1 + sqrt(2)) .* V_edf.pb ./ (1 + (4 + V_edf.pb.^4).^0.25);
    W_edf.pb               = (V_edf.pb.^2 - U_edf.pb.^2).^0.5;
    w_edf.pb               = r_edf .* V_edf.pb .* besselk(1, W_edf.pb)...  % модова€ площадь (встречна€ накачка)
    .* besselj(0, U_edf.pb) ./ (U_edf.pb .* besselk(0, W_edf.pb));
psi.pb(1: N.pb,1: n_r)   = (besselj(0, (repmat(U_edf.pb',1,n_r)...
    .* repmat(r, N.pb,1)) / r_edf)).^2;
psi.npb            = psi.pb ./ pi ./ w_edf.pb'.^2;                     % нормированна€ функци€ распределени€ (встречна€ накачка)
end
% дл€  r>a
% r_more_a = r_edf: dr : 10^(-5);
% n_r_ma   = size(r_more_a,2);                                          % размер массива
% 
% psi.S(1: N.S,n_r_la+1: n_r_la+n_r_ma)     = (besselk(0,(repmat(W_edf.S',1,n_r_ma) .*...
%     repmat(r_more_a, N.S,1)) / r_edf)).^2 .* besselj(0, repmat(U_edf.S',1,n_r_ma)).^2 ./...
%     besselk(0, repmat(W_edf.S',1,n_r_ma)).^2;  
% psi.PF(1: N.PF,n_r_la+1: n_r_la+n_r_ma)   = (besselk(0,(repmat(W_edf.PF',1,n_r_ma)...
%     .* repmat(r_more_a, N.PF,1)) / r_edf)).^2 .* besselj(0, repmat(U_edf.PF',1,n_r_ma)).^2 ./...
%     besselk(0, repmat(W_edf.PF',1,n_r_ma)).^2;
% psi.PB(1: N.PB,n_r_la+1: n_r_la+n_r_ma)   = (besselk(0,(repmat(W_edf.PB',1,n_r_ma)...
%     .* repmat(r_more_a, N.PB,1)) / r_edf)).^2 .* besselj(0, repmat(U_edf.PB',1,n_r_ma)).^2 ./...
%     besselk(0, repmat(W_edf.PB',1,n_r_ma)).^2;
% psi.ASE(1: N.ASE,n_r_la+1: n_r_la+n_r_ma) = (besselk(0,(repmat(W_edf.ASE',1,n_r_ma)...
%     .* repmat(r_more_a, N.ASE,1)) / r_edf)).^2 .* besselj(0, repmat(U_edf.ASE',1,n_r_ma)).^2 ./...
%     besselk(0, repmat(W_edf.ASE',1,n_r_ma)).^2;
% psi.NS              = psi.S ./ pi ./ w_edf.S'.^2;    % нормированна€ функци€ распределени€ (сигнал)
% psi.NPF             = psi.PF ./ pi ./ w_edf.PF'.^2;  % нормированна€ функци€ распределени€ (попутна€ накачка)
% psi.NPB             = psi.PB ./ pi ./ w_edf.PB'.^2;  % нормированна€ функци€ распределени€ (встречна€ накачка)
% psi.NASE            = psi.ASE ./ pi ./ w_edf.ASE'.^2;% нормированна€ функци€ распределени€ (ASE)

end