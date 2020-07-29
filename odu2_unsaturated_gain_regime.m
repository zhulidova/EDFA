function f = odu2_unsaturated_gain_regime(z, P, wl, sigma, N,w_edf, n_sum, ph_const, P_in)
%% упрощенная система ОДУ (режим слабого сигнала)
% расчет интегралов при замене волновых функций интенсивности экспонентами (вывод -
%отчеты по модели EDFA/Жулидова. Отчет от 02.06)

% Переменные
% P_sat   - структура мощностей насыщения для сигнала (P_sat.S), накачки (P_sat.PF и P_sat.PB) и ASE (P_sat.ASE)
% p       - структура приведенных мощностей сигнала (p.S), накачки (p.PF и p.PB) и ASE (p.ASEF и p.ASEB)
% overlap - структура интегралов перекрытия: усиление (overlap.ES, overlap.EASE), поглощение (overlap.AS, overlap.AASE)
%overlap.sum_X - переменные для расчета интегралов перекрытия
% g       - структура дифференциальных коэффициентов усиления (g.ES и g.EASE) и поглощения (g.AS и g.AASE)
% alpha   - структура коэффицентов истощения попутной (alpha.PF) и встречной (alpha.PB) накачек
% f       - значение производных в скоростных уравнениях

P_sat.S             = (ph_const.h * ph_const.c * 10^34 * pi .* w_edf.S.^2)...   % мощность насыщения (сигнал)
    ./ (wl.S .* (sigma.AS + sigma.ES) * ph_const.tau);
P_sat.PF            = (ph_const.h * ph_const.c * 10^34 * pi .* w_edf.PF.^2)...  % мощность насыщения (попутная накачка)
    ./ (wl.PF .* (sigma.APF + sigma.EPF) * ph_const.tau);                                   
P_sat.PB            = (ph_const.h * ph_const.c * 10^34 * pi .* w_edf.PB.^2)...  % мощность насыщения (встречная накачка)
    ./ (wl.PB .* (sigma.APB + sigma.EPB) * ph_const.tau);                                   
P_sat.ASE           = (ph_const.h * ph_const.c * 10^34 * pi .* w_edf.ASE.^2)... % мощность насыщения (ASE)
    ./ (wl.ASE .* (sigma.AASE + sigma.EASE) * ph_const.tau);

p.S                 = P(1: N.S, 1) ./ P_sat.S';                                 % приведенная мощность сигнала
p.ASEF              = P(N.S+1:N.S+N.ASE, 1) ./ P_sat.ASE';                      % приведенная мощность попутного ASE
p.PF                = P(N.S+N.ASE+1:N.S+N.ASE+N.PF,1) ./ P_sat.PF';             % приведенная мощность попутной накачки
overlap.sum_ES      = sum(p.S' .* sigma.ES ./ (sigma.AS + sigma.ES)) +...       % слагаемое с усилением для сигнала и ASE
    sum(p.ASEF' .* sigma.EASE ./ (sigma.AASE + sigma.EASE));
overlap.sum_all     = sum(p.S) + sum(p.ASEF) + sum(p.PF);                       % сумма всех приведенных мощностей
overlap.sum_EP      = sum(p.PF' .* sigma.EPF ./ (sigma.APF + sigma.EPF));      % слагаемое с усилением для накачки
overlap.sum_AP      = sum(p.PF' .* sigma.APF ./ (sigma.APF + sigma.EPF));       % слагаемое с поглощением для накачки

if isempty(P_in.PB) == 0
    p.ASEB          = P(N.S+N.ASE+N.PF+1:N.S+N.ASE+N.PF+N.ASE) ./ P_sat.ASE';   % приведенная мощность попутного ASE
    p.PB            = P(N.S+N.ASE+N.PF+N.ASE+1:N.S+N.ASE+N.PF+N.ASE+N.PB) ...   % приведенная мощность встречной накачки
        ./ P_sat.PB';
    overlap.sum_ES  = overlap.sum_ES + sum(p.ASEB' .* sigma.EASE...             % слагаемое с усилением для сигнала и ASE
        ./ (sigma.AASE + sigma.EASE));
    overlap.sum_EP  = overlap.sum_EP + sum(p.PB' .* sigma.EPB ./ ...            % слагаемое с усилением для накачки
        (sigma.APB + sigma.EPB));
    overlap.sum_AP  = overlap.sum_AP + sum(p.PB' .* sigma.APB ./ ...            % слагаемое с поглощением для накачки
        (sigma.APB + sigma.EPB));
    overlap.sum_all = overlap.sum_all + sum(p.ASEB) + sum(p.PB);                % сумма всех приведенных мощностей
end

overlap.ES      = overlap.sum_AP ./ overlap.sum_all.^2 .* (overlap.sum_all...  % интеграл перекрытия (усиление) для сигнала
    .* (1 - exp(-(10^(-5) ./ w_edf.S).^2)));
overlap.AS      = (overlap.sum_ES + overlap.sum_EP) ./ overlap.sum_all.^2 .*...% интеграл перекрытия (поглощение) для сигнала
    (overlap.sum_all .* (1 - exp(-(10^(-5) ./ w_edf.S).^2)));
overlap.EASE    = overlap.sum_AP ./ overlap.sum_all.^2 .* (overlap.sum_all .*...% интеграл перекрытия (усиление) для ASE
    (1 - exp(-(10^(-5) ./ w_edf.ASE).^2)));
overlap.AASE    = (overlap.sum_ES + overlap.sum_EP) ./ overlap.sum_all.^2 .*... % интеграл перекрытия (поглощение) для ASE
    (overlap.sum_all .* (1 - exp(-(10^(-5) ./ w_edf.ASE).^2)));

g.ES            = n_sum .* sigma.ES .* overlap.ES;                              % диф. коэффициент усиления (сигнал)
g.AS            = n_sum .* sigma.AS .* overlap.AS;                              % диф. коэффициент поглощения (сигнал)
g.EASE          = n_sum .* sigma.EASE .* overlap.EASE;                          % диф. коэффициент усиления (ASE)
g.AASE          = n_sum .* sigma.AASE .* overlap.AASE;                          % диф. коэффициент поглощения (ASE)

alpha.PF        = (1 - exp(-(10^(-5) ./ w_edf.PF).^2)) * n_sum .* sigma.APF;    % коэффициент истощения попутной накачки
alpha.PB        = (1 - exp(-(10^(-5) ./ w_edf.PB).^2)) * n_sum .* sigma.APB;    % коэффицент истощения встречной накачки

% итоговые значения дифференциалов

%уравнения для сигнала (dP_s/dz)
f(1:N.S,1)  = (g.ES - g.AS)' .* P(1: N.S, 1);

% уравнения для попутного ASE (dP_ase/dz)
f(N.S+1:N.S+N.ASE,1)  = (g.EASE - g.AASE)' .* P(N.S+1:N.S+N.ASE, 1) + 2 * ph_const.h * ph_const.c^2 * 0.1 *...
    10^18 ./ wl.ASE'.^3 .* g.EASE';

% уравнения для попутной накачки (dP_p/dz)
f(N.S+N.ASE+1:N.S+N.ASE+N.PF,1)  = -alpha.PF' ./ (1 + sum(p.PF)) .* P(N.S+N.ASE+1:N.S+N.ASE+N.PF,1);

% допольнительные уравнения при встречной накачке
if isempty(P_in.PB) == 0
    % уравнения для встречного ASE (dP_ase/dz)
    f(N.S+N.ASE+N.PF+1:N.S+N.ASE+N.PF+N.ASE,1)  = (-1) .* (g.EASE - g.AASE)' .*...
        P(N.S+N.ASE+N.PF+1: N.S+N.ASE+N.PF+N.ASE, 1) - 2 * ph_const.h * ph_const.c^2 * 0.1 * 10^18 ./ wl.ASE'.^3 ...
        .* g.EASE';
    
    % уравнение для встречной накачки (dP_p/dz)
    f(N.S+N.ASE+N.PF+N.ASE+1:N.S+N.ASE+N.PF+N.ASE+N.PB,1)  = alpha.PB' ./ (1 + sum(p.PB)) .*...
        P(N.S+N.ASE+N.PF+N.ASE+1: N.S+N.ASE+N.PF+N.ASE+N.PB, 1);   
end
end