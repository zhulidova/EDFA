%% расчет населенностей верхнего и нижнего энергетических уровней
function [n1, n2] = population(P, wl, sigma, psi, N, ph_const, P_in, r_edf)
% Параметры 
% P_sat   - структура значений мощност насыщения сигнала(P_sat.S), накачки (P_sat.PF и P_sat.PB) и ASE (P_sat.ASE)
% W12     - вероятность вынужденного перехода с поглощением излучения для сигнала
% W21     - вероятность вынужденного излучательного перехода для сигнала
% R12     - вероятность вынужденного перехода с поглощением излучения для накачки
% R21     - вероятность вынужденного излучательного перехода для накачки
% n1 и n2 - населенности верхнего и нижнего уровня

P_sat.S   = (ph_const.h * ph_const.c * 10^34) ./ (wl.S .* (sigma.AS + sigma.ES));                  % мощность насыщения для сигнала
P_sat.ASE = (ph_const.h * ph_const.c * 10^34) ./ (wl.ASE .* (sigma.AASE + sigma.EASE));            % мощность насыщения для ASE
P_sat.PF  = (ph_const.h * ph_const.c * 10^34) ./ (wl.PF .* (sigma.APF + sigma.EPF));               % мощность насыщения для попутной накачки
% вероятность вынужденного перехода с поглощением излучения для сигнала и ASE
W12       = sum(psi.NS .* repmat(P(1: N.S, 1),1,size(psi.NS,2)) .* repmat(sigma.AS',1,size(psi.NS,2))...
    ./ ((repmat(sigma.AS',1,size(psi.NS,2)) + repmat(sigma.ES',1,size(psi.NS,2))) .* repmat(P_sat.S',1,size(psi.NS,2))))+...
    sum(psi.NASE .* repmat(P(N.S+1: N.S+N.ASE, 1),1,size(psi.NASE,2)) .* repmat(sigma.AASE',1,size(psi.NASE,2))...
    ./ ((repmat(sigma.AASE',1,size(psi.NASE,2)) + repmat(sigma.EASE',1,size(psi.NASE,2))) .* repmat(P_sat.ASE',1,size(psi.NASE,2))));
% вероятность вынужденного излучательного перехода для сигнала и ASE 
W21       = sum(psi.NS .* repmat(P(1: N.S, 1),1,size(psi.NS,2)) .* repmat(sigma.ES',1,size(psi.NS,2))...
    ./ ((repmat(sigma.AS',1,size(psi.NS,2)) + repmat(sigma.ES',1,size(psi.NS,2))) .* repmat(P_sat.S',1,size(psi.NS,2))))+...
    sum(psi.NASE .* repmat(P(N.S+1: N.S+N.ASE, 1),1,size(psi.NASE,2)) .* repmat(sigma.EASE',1,size(psi.NASE,2))...
    ./ ((repmat(sigma.AASE',1,size(psi.NASE,2)) + repmat(sigma.EASE',1,size(psi.NASE,2))) .* repmat(P_sat.ASE',1,size(psi.NASE,2))));   
% вероятность вынужденного перехода с поглощением излучения для накачки
R12    = sum(repmat(P(N.S+N.ASE+1: N.S+N.ASE+N.PF, 1),1,size(psi.NPF,2)) .* psi.NPF .* repmat(sigma.APF',1,size(psi.NPF,2))...
    ./ ((repmat(sigma.APF',1,size(psi.NPF,2)) + repmat(sigma.EPF',1,size(psi.NPF,2))) .* repmat(P_sat.PF',1,size(psi.NPF,2))));
% вероятность вынужденного излучательного перехода для накачки
R21    = sum(repmat(P(N.S+N.ASE+1: N.S+N.ASE+N.PF, 1),1,size(psi.NPF,2)) .* psi.NPF .* repmat(sigma.EPF',1,size(psi.NPF,2))...
    ./ ((repmat(sigma.APF',1,size(psi.NPF,2)) + repmat(sigma.EPF',1,size(psi.NPF,2))) .* repmat(P_sat.PF',1,size(psi.NPF,2))));

% допольнительные уравнения при встречной накачке
if isempty(P_in.PB) == 0 
    P_sat.PB = (ph_const.h * ph_const.c * 10^34) ./ (wl.PB .* (sigma.APB + sigma.EPB));
    % вероятность вынужденного излучательного перехода для сигнала и ASE
    W21      = W21 + sum(psi.NASE .* repmat(P(N.S+N.ASE+N.PF+1: N.S+N.ASE+N.ASE+N.PF, 1),1,size(psi.NASE,2)) .*...
        repmat(sigma.EASE',1,size(psi.NASE,2)) ./ ((repmat(sigma.AASE',1,size(psi.NASE,2)) +...
        repmat(sigma.EASE',1,size(psi.NASE,2))) .* repmat(P_sat.ASE',1,size(psi.NASE,2))));
    % вероятность вынужденного перехода с поглощением излучения для накачки
    R12      = R12 + sum(repmat(P(N.S+N.ASE+N.ASE+N.PF+1: N.S+N.ASE+N.ASE+N.PF+N.PB, 1),1,size(psi.NPB,2)) .*...
        psi.NPB .* repmat(sigma.APB',1,size(psi.NPB,2)) ./ ((repmat(sigma.APB',1,size(psi.NPB,2)) +...
        repmat(sigma.EPB',1,size(psi.NPB,2))) .* repmat(P_sat.PB',1,size(psi.NPB,2))));
    % вероятность вынужденного излучательного перехода для накачки
    R21      = R21 + sum(repmat(P(N.S+N.ASE+N.ASE+N.PF+1: N.S+N.ASE+N.ASE+N.PF+N.PB, 1),1,size(psi.NPB,2)) .*...
        psi.NPB .* repmat(sigma.EPB',1,size(psi.NPB,2)) ./ ((repmat(sigma.APB',1,size(psi.NPB,2)) +...
        repmat(sigma.EPB',1,size(psi.NPB,2))) .* repmat(P_sat.PB',1,size(psi.NPB,2))));
end

% населенности верхнего и нижнего уровня    
n2           = (W12 + R12) ./ (W12 + R12 + W21 + 1 / ph_const.tau + R21);                           % населенность верхнего уровня
n1           = (1 / ph_const.tau + W21 + R21) ./ (W12 + R12 + W21 + 1 / ph_const.tau + R21);        % населенность нижнего уровня

% концетрация ионов эрбия при r>r_edf равна нулю
r_length_max = length(r_edf: 10^(-7): 10^(-5));                                                     % длина массива r>r_edf
r_length_min = length(0: 10^(-7): r_edf);                                                           % длина массива r<r_edf
n1(1, r_length_min: r_length_max) = zeros(1, r_length_max - r_length_min + 1);
n2(1, r_length_min: r_length_max) = zeros(1, r_length_max - r_length_min + 1);
end