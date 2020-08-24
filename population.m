%% расчет населенностей верхнего и нижнего энергетических уровней
function [n1, n2] = population(P, Lambda, sigma, psi, N, ph_const, Pin)
% Параметры 
% P_sat   - структура значений мощност насыщения сигнала(P_sat.S), накачки (P_sat.PF и P_sat.PB) и ASE (P_sat.ASE)
% W12     - вероятность вынужденного перехода с поглощением излучения для сигнала
% W21     - вероятность вынужденного излучательного перехода для сигнала
% R12     - вероятность вынужденного перехода с поглощением излучения для накачки
% R21     - вероятность вынужденного излучательного перехода для накачки
% n1 и n2 - населенности верхнего и нижнего уровня

P_sat.s   = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.s .* (sigma.as + sigma.es));                  % мощность насыщения для сигнала
P_sat.ase = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.ase .* (sigma.aase + sigma.ease));            % мощность насыщения для ASE
P_sat.pf  = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.pf .* (sigma.apf + sigma.epf));               % мощность насыщения для попутной накачки
% вероятность вынужденного перехода с поглощением излучения для сигнала и ASE
W12       = sum(psi.ns .* repmat(P(1: N.s, 1),1,size(psi.ns,2)) .* repmat(sigma.as',1,size(psi.ns,2))...
    ./ ((repmat(sigma.as',1,size(psi.ns,2)) + repmat(sigma.es',1,size(psi.ns,2))) .* repmat(P_sat.s',1,size(psi.ns,2))))+...
    sum(psi.nase .* repmat(P(N.s+1: N.s+N.ase, 1),1,size(psi.nase,2)) .* repmat(sigma.aase',1,size(psi.nase,2))...
    ./ ((repmat(sigma.aase',1,size(psi.nase,2)) + repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))));
% вероятность вынужденного излучательного перехода для сигнала и ASE 
W21       = sum(psi.ns .* repmat(P(1: N.s, 1),1,size(psi.ns,2)) .* repmat(sigma.es',1,size(psi.ns,2))...
    ./ ((repmat(sigma.as',1,size(psi.ns,2)) + repmat(sigma.es',1,size(psi.ns,2))) .* repmat(P_sat.s',1,size(psi.ns,2))))+...
    sum(psi.nase .* repmat(P(N.s+1: N.s+N.ase, 1),1,size(psi.nase,2)) .* repmat(sigma.ease',1,size(psi.nase,2))...
    ./ ((repmat(sigma.aase',1,size(psi.nase,2)) + repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))));   
% вероятность вынужденного перехода с поглощением излучения для накачки
R12    = sum(repmat(P(N.s+N.ase+1: N.s+N.ase+N.pf, 1),1,size(psi.npf,2)) .* psi.npf .* repmat(sigma.apf',1,size(psi.npf,2))...
    ./ ((repmat(sigma.apf',1,size(psi.npf,2)) + repmat(sigma.epf',1,size(psi.npf,2))) .* repmat(P_sat.pf',1,size(psi.npf,2))));
% вероятность вынужденного излучательного перехода для накачки
R21    = sum(repmat(P(N.s+N.ase+1: N.s+N.ase+N.pf, 1),1,size(psi.npf,2)) .* psi.npf .* repmat(sigma.epf',1,size(psi.npf,2))...
    ./ ((repmat(sigma.apf',1,size(psi.npf,2)) + repmat(sigma.epf',1,size(psi.npf,2))) .* repmat(P_sat.pf',1,size(psi.npf,2))));

% допольнительные уравнения при встречной накачке
if isempty(Pin.pb) == 0 
    P_sat.pb = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.pb .* (sigma.apb + sigma.epb));
    % вероятность вынужденного излучательного перехода для сигнала и ASE
    W21      = W21 + sum(psi.nase .* repmat(P(N.s+N.ase+N.pf+1: N.s+2*N.ase+N.pf, 1),1,size(psi.nase,2)) .*...
        repmat(sigma.ease',1,size(psi.nase,2)) ./ ((repmat(sigma.aase',1,size(psi.nase,2)) +...
        repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))));
    % вероятность вынужденного перехода с поглощением излучения для накачки
    R12      = R12 + sum(repmat(P(N.s+2*N.ase+N.pf+1: N.s+2*N.ase+N.pf+N.pb, 1),1,size(psi.npb,2)) .*...
        psi.npb .* repmat(sigma.apb',1,size(psi.npb,2)) ./ ((repmat(sigma.apb',1,size(psi.npb,2)) +...
        repmat(sigma.epb',1,size(psi.npb,2))) .* repmat(P_sat.pb',1,size(psi.npb,2))));
    % вероятность вынужденного излучательного перехода для накачки
    R21      = R21 + sum(repmat(P(N.s+2*N.ase+N.pf+1: N.s+2*N.ase+N.pf+N.pb, 1),1,size(psi.npb,2)) .*...
        psi.npb .* repmat(sigma.epb',1,size(psi.npb,2)) ./ ((repmat(sigma.apb',1,size(psi.npb,2)) +...
        repmat(sigma.epb',1,size(psi.npb,2))) .* repmat(P_sat.pb',1,size(psi.npb,2))));
end

% населенности верхнего и нижнего уровня    
n2           = (W12 + R12) ./ (W12 + R12 + W21 + 1 / ph_const.tau + R21);                           % населенность верхнего уровня
n1           = (1 / ph_const.tau + W21 + R21) ./ (W12 + R12 + W21 + 1 / ph_const.tau + R21);        % населенность нижнего уровня

end