%% Функция расчета входной мощности ASE
%(Desurvire, appendix R)
% P_sat_P   - массив значений мощности насыщения для накачки
% P_sat_ASE - массив значений мощности насыщения для спектра ASE
% n_sp      - массив коэффициентов спонтанной люминесценции, имеет размер массива ASE
% p_ase     - мощность ASE в Вт
% P_ase     - мощность ASE в дБм

function P_ase = ase_in(Lambda, ph_const, sigma, PinP, w_edf) 

P_sat_P   = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.pf.^2)...                  % мощность насыщения для накачки
    ./ (Lambda.pf .* (sigma.apf + sigma.epf) * ph_const.tau);
P_sat_ASE = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.ase.^2)...                 % мощность насыщения для ASE
    ./ (Lambda.ase .* (sigma.aase + sigma.ease) * ph_const.tau);

n_sp      = 1 ./ (1 - sum(sigma.epf ./ sigma.apf) .* sigma.aase  ./ sigma.ease...     % коэффициент спонтанной люминесценции
    - sum((1 + sigma.epf./ sigma.apf) .* P_sat_P ./ PinP) .* sigma.aase  ./ sigma.ease);

if n_sp(1,1) < 0
    n_sp      = 1 ./ (1 - sum(sigma.epf ./ sigma.apf) .* sigma.aase  ./ sigma.ease...     % коэффициент спонтанной люминесценции
    - sum((1 + sigma.epf./ sigma.apf) .* P_sat_P ./ (3.5.*PinP)) .* sigma.aase  ./ sigma.ease);
end
    
P_ase     = 4 * n_sp .* ph_const.c * 0.1 * 10^(-9) ./ Lambda.ase.^2 .* P_sat_ASE * 10^(-16); % мощность ASE в Вт

end