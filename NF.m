%% Расчет спектра шум-фактора
% расчет из Giles, Desurvire, 1991, Propagation of signal and noise in
% concatenated erbium-doped fiber optical amplifiers

% Формула шум-фактора (Desurvire):
% (1 + 2 .* n_sp .* (undb(G) - 1)) ./ undb(G)            - дробовой шум
% M .* (n_sp .* (undb(G) - 1)).^2 ./(undb(G)).^2/n_mean  - спонтанно-спонтанный шум биений
% M .* n_sp .* (undb(G) - 1)./(undb(G)).^2/n_mean        - сигнал-спонтанный шум биений

function [OSNR, NF] = nf(P_out_ase,G, Lambda, ph_const, P_out_s)

P_ase  = interp1(Lambda.ase,P_out_ase',Lambda.s);          % вычисление мощности ASE на длине волны сигнала

% M      = 1;                                              % число мод
% n_mean = 2;                                              % среднее число фотонов при z = 0, распределение Пуассона

n_sp   = P_ase' ./ (2 * (undb(G) - 1) * ph_const.h...      % коэффициент спонтанной люминесценции
    * ph_const.c^2 * 0.1 * 10^(-9) ./ Lambda.s'.^3);   
f      = (1 + 2 .* n_sp .* (undb(G) - 1)) ./ undb(G);      % шум-фактор в линейных единицах
NF     = db(f);                                            % шум-фактор в дБ
OSNR   = db(P_out_s ./ P_ase');
end