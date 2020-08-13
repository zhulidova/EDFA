%% Расчет спектра шум-фактора
% расчет из Giles, Desurvire, 1991, Propagation of signal and noise in
% concatenated erbium-doped fiber optical amplifiers

% Формула шум-фактора (Desurvire):
% (1 + 2 .* n_sp .* (undb(G) - 1)) ./ undb(G)            - дробовой шум
% M .* (n_sp .* (undb(G) - 1)).^2 ./(undb(G)).^2/n_mean  - спонтанно-спонтанный шум биений
% M .* n_sp .* (undb(G) - 1)./(undb(G)).^2/n_mean        - сигнал-спонтанный шум биений

function nf = NF(P_out_ase,G, wl, ph_const)

P_ase  = interp1(wl.ASE,P_out_ase,wl.S);                   % вычисление мощности ASE на длине волны сигнала

M      = 1;                                                % число мод
n_mean = 2;                                                % среднее число фотонов при z = 0, распределение Пуассона

n_sp   = P_ase ./ (2 * (undb(G) - 1) * ph_const.h...       % коэффициент спонтанной люминесценции
    * ph_const.c^2 * 0.1 * 10^18 ./ wl.S.^3);   
f      = (1 + 2 .* n_sp .* (undb(G) - 1)) ./ undb(G);      % шум-фактор в линейных единицах
nf     = db(f);                                            % шум-фактор в дБ
end