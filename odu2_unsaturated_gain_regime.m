function f = odu2_unsaturated_gain_regime(z, P, Lambda, sigma, N,w_edf, n_sum, ph_const, P_in)
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

P_sat.s             = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.s.^2)...   % мощность насыщения (сигнал)
    ./ (Lambda.s .* (sigma.as + sigma.es) * ph_const.tau);
P_sat.pf            = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.pf.^2)...  % мощность насыщения (попутная накачка)
    ./ (Lambda.pf .* (sigma.apf + sigma.epf) * ph_const.tau);                                   
P_sat.pb            = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.pb.^2)...  % мощность насыщения (встречная накачка)
    ./ (Lambda.pb .* (sigma.apb + sigma.epb) * ph_const.tau);                                   
P_sat.ase           = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.ase.^2)... % мощность насыщения (ASE)
    ./ (Lambda.ase .* (sigma.aase + sigma.ease) * ph_const.tau);

p.s                 = P(1: N.s, 1) ./ P_sat.s';                                 % приведенная мощность сигнала
p.asef              = P(N.s+1:N.s+N.ase, 1) ./ P_sat.ase';                      % приведенная мощность попутного ASE
p.pf                = P(N.s+N.ase+1:N.s+N.ase+N.pf,1) ./ P_sat.pf';             % приведенная мощность попутной накачки
overlap.sum_es      = sum(p.s' .* sigma.es ./ (sigma.as + sigma.es)) +...       % слагаемое с усилением для сигнала и ASE
    sum(p.asef' .* sigma.ease ./ (sigma.aase + sigma.ease));
overlap.sum_all     = sum(p.s) + sum(p.asef) + sum(p.pf);                       % сумма всех приведенных мощностей
overlap.sum_ep      = sum(p.pf' .* sigma.epf ./ (sigma.apf + sigma.epf));       % слагаемое с усилением для накачки
overlap.sum_ap      = sum(p.pf' .* sigma.apf ./ (sigma.apf + sigma.epf));       % слагаемое с поглощением для накачки

if isempty(P_in.pb) == 0
    p.aseb          = P(N.s+N.ase+N.pf+1:N.s+2*N.ase+N.pf) ./ P_sat.ase';   % приведенная мощность попутного ASE
    p.pb            = P(N.s+2*N.ase+N.pf+1:N.s+2*N.ase+N.pf+N.pb) ...   % приведенная мощность встречной накачки
        ./ P_sat.pb';
    overlap.sum_es  = overlap.sum_es + sum(p.aseb' .* sigma.ease...             % слагаемое с усилением для сигнала и ASE
        ./ (sigma.aase + sigma.ease));
    overlap.sum_ep  = overlap.sum_ep + sum(p.pb' .* sigma.epb ./ ...            % слагаемое с усилением для накачки
        (sigma.apb + sigma.epb));
    overlap.sum_ap  = overlap.sum_ap + sum(p.pb' .* sigma.apb ./ ...            % слагаемое с поглощением для накачки
        (sigma.apb + sigma.epb));
    overlap.sum_all = overlap.sum_all + sum(p.aseb) + sum(p.pb);                % сумма всех приведенных мощностей
end

overlap.es      = overlap.sum_ap ./ overlap.sum_all.^2 .* (overlap.sum_all...   % интеграл перекрытия (усиление) для сигнала
    .* (1 - exp(-(10^(-5) ./ w_edf.s).^2)));
overlap.as      = (overlap.sum_es + overlap.sum_ep) ./ overlap.sum_all.^2 .*... % интеграл перекрытия (поглощение) для сигнала
    (overlap.sum_all .* (1 - exp(-(10^(-5) ./ w_edf.s).^2)));
overlap.ease    = overlap.sum_ap ./ overlap.sum_all.^2 .* (overlap.sum_all .*...% интеграл перекрытия (усиление) для ASE
    (1 - exp(-(10^(-5) ./ w_edf.ase).^2)));
overlap.aase    = (overlap.sum_es + overlap.sum_ep) ./ overlap.sum_all.^2 .*... % интеграл перекрытия (поглощение) для ASE
    (overlap.sum_all .* (1 - exp(-(10^(-5) ./ w_edf.ase).^2)));

g.es            = n_sum .* sigma.es .* overlap.es;                              % диф. коэффициент усиления (сигнал)
g.as            = n_sum .* sigma.as .* overlap.as;                              % диф. коэффициент поглощения (сигнал)
g.ease          = n_sum .* sigma.ease .* overlap.ease;                          % диф. коэффициент усиления (ASE)
g.aase          = n_sum .* sigma.aase .* overlap.aase;                          % диф. коэффициент поглощения (ASE)

alpha.pf        = (1 - exp(-(10^(-5) ./ w_edf.pf).^2)) * n_sum .* sigma.apf;    % коэффициент истощения попутной накачки
alpha.pb        = (1 - exp(-(10^(-5) ./ w_edf.pb).^2)) * n_sum .* sigma.apb;    % коэффицент истощения встречной накачки

% итоговые значения дифференциалов

%уравнения для сигнала (dP_s/dz)
f(1:N.s,1)  = (g.es - g.as)' .* P(1: N.s, 1);

% уравнения для попутного ASE (dP_ase/dz)
f(N.s+1:N.s+N.ase,1)  = (g.ease - g.aase)' .* P(N.s+1:N.s+N.ase, 1) + 2 * ph_const.h * ph_const.c^2 * 0.1 *...
    10^(-9) ./ Lambda.ase'.^3 .* g.ease';

% уравнения для попутной накачки (dP_p/dz)
f(N.s+N.ase+1:N.s+N.ase+N.pf,1)  = -alpha.pf' ./ (1 + sum(p.pf)) .* P(N.s+N.ase+1:N.s+N.ase+N.pf,1);

% допольнительные уравнения при встречной накачке
if isempty(P_in.pb) == 0
    % уравнения для встречного ASE (dP_ase/dz)
    f(N.s+N.ase+N.pf+1:N.s+2*N.ase+N.pf,1)  = (-1) .* (g.ease - g.aase)' .*...
        P(N.s+N.ase+N.pf+1: N.s+2*N.ase+N.pf, 1) - 2 * ph_const.h * ph_const.c^2 * 0.1 * 10^(-9) ./ Lambda.ase'.^3 ...
        .* g.ease';
    
    % уравнение для встречной накачки (dP_p/dz)
    f(N.s+2*N.ase+N.pf+1:N.s+2*N.ase+N.pf+N.pb,1)  = alpha.pb' ./ (1 + sum(p.pb)) .*...
        P(N.s+2*N.ase+N.pf+1: N.s+2*N.ase+N.pf+N.pb, 1);   
end
end