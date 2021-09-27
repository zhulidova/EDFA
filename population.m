%% расчет населенностей верхнего и нижнего энергетических уровней
function n = population(P, Lambda, sigma, psi, N, ph_const,n_sum)
% Параметры 
% P_sat   - структура значений мощност насыщения сигнала(P_sat.s), накачки (P_sat.pf и P_sat.pb) и ASE (P_sat.ase)
% W12     - скорость вынужденного перехода с поглощением излучения для сигнала
% W21     - скорость вынужденного излучательного перехода для сигнала
% R12     - скорость вынужденного перехода с поглощением излучения для накачки
% R21     - скорость вынужденного излучательного перехода для накачки
% R24     - скорость вынужденного перехода с поглощением с уровня I(13/2)
% на I(9/2), длина волны 1680 и люминесценции с (11/2) на F(7/2), длина волны 980
% n1      - массив населенностей уровня 1 : I(15/2), при разном положении в сердцевине
% n2      - массив населенностей уровня 2 : I(13/2), при разном положении в сердцевине
% n3      - массив населенностей уровня 3 : I(11/2), при разном положении в сердцевине
% n4      - массив населенностей уровня 4 : I(19/2), при разном положении в сердцевине

P_sat.s   = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.s .* (sigma.as + sigma.es));       % мощность насыщения для сигнала
P_sat.ase = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.ase .* (sigma.aase + sigma.ease)); % мощность насыщения для ASE
P_sat.pf  = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.pf .* (sigma.apf + sigma.epf));    % мощность насыщения для попутной накачки

R12 = 0;
R21 = 0;
R13 = 0;
R31 = 0;
R24 = 0;
n.first  = zeros(1,size(psi.ns,2));
n.second = zeros(1,size(psi.ns,2));
n.third  = zeros(1,size(psi.ns,2));
n.fourth = zeros(1,size(psi.ns,2));

%% Рассчет скоростей переходов
% скорость вынужденного перехода с поглощением излучения для сигнала и ASE
W12       = sum(psi.ns .* repmat(P(1: N.s, 1),1,size(psi.ns,2)) .* repmat(sigma.as',1,size(psi.ns,2))...
    ./ ((repmat(sigma.as',1,size(psi.ns,2)) + repmat(sigma.es',1,size(psi.ns,2))) .* repmat(P_sat.s',1,size(psi.ns,2))))+...
    sum(psi.nase .* repmat(P(N.s+1: N.s+N.ase, 1),1,size(psi.nase,2)) .* repmat(sigma.aase',1,size(psi.nase,2))...
    ./ ((repmat(sigma.aase',1,size(psi.nase,2)) + repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))),1);

% скорость вынужденного излучательного перехода для сигнала и ASE 
W21       = sum(psi.ns .* repmat(P(1: N.s, 1),1,size(psi.ns,2)) .* repmat(sigma.es',1,size(psi.ns,2))...
    ./ ((repmat(sigma.as',1,size(psi.ns,2)) + repmat(sigma.es',1,size(psi.ns,2))) .* repmat(P_sat.s',1,size(psi.ns,2))))+...
    sum(psi.nase .* repmat(P(N.s+1: N.s+N.ase, 1),1,size(psi.nase,2)) .* repmat(sigma.ease',1,size(psi.nase,2))...
    ./ ((repmat(sigma.aase',1,size(psi.nase,2)) + repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))),1); 

R24  = sum(psi.ns .* repmat(P(1: N.s, 1),1,size(psi.ns,2)) .* repmat(sigma.esa_eta.*sigma.as',1,size(psi.ns,2)) ./...
    (ph_const.h * ph_const.c * 10^25).* repmat(Lambda.s',1,size(psi.ns,2)),1);

R24  = R24 + sum(psi.nase .* repmat(P(N.s+1: N.s+N.ase, 1),1,size(psi.nase,2)) .* repmat(sigma.esa_eta.*sigma.aase',1,size(psi.nase,2))...
    ./ (ph_const.h * ph_const.c * 10^25).*repmat(Lambda.ase',1,size(psi.nase,2)),1);


for i = 1:length(Lambda.pf)
if Lambda.pf(i) > 1450* 10^(-9) % для накачки с длиной волны 1480 нм
    
    % скорость вынужденного перехода с поглощением излучения для накачки
    R12  = R12 + psi.npf(i,:) .* P(N.s+N.ase+i, 1) .* sigma.apf(i) ./ (sigma.apf(i) + sigma.epf(i)) ./ P_sat.pf(i);
    
    % скорость вынужденного излучательного перехода для накачки
    R21  = R21 + psi.npf(i,:) .* P(N.s+N.ase+i, 1) .* sigma.epf(i) ./ (sigma.apf(i) + sigma.epf(i)) ./ P_sat.pf(i);
    
    % учет эффекта ESA для накачки на длине волны 1480 нм (для 980 нм должна быть большая мощность накачки)
    R24  = R24 + psi.npf(i,:) .* P(N.s+N.ase+i, 1) .* sigma.esa_eta*sigma.apf(i) ./ (ph_const.h * ph_const.c * 10^25).*Lambda.pf(i);


else                            % для накачки с длиной волны 980 нм
    
    % скорость вынужденного перехода с поглощением излучения для накачки
    R13  = R13 + psi.npf(i,:) .* P(N.s+N.ase+i, 1) .* sigma.apf(i) ./ (sigma.apf(i) + sigma.epf(i)) ./ P_sat.pf(i);
    
    % скорость вынужденного излучательного перехода для накачки
    R31  = R31 + psi.npf(i,:) .* P(N.s+N.ase+i, 1) .* sigma.epf(i) ./ (sigma.apf(i) + sigma.epf(i)) ./ P_sat.pf(i);
end
end
% допольнительные уравнения при встречной накачке
if isempty(Lambda.pb) == 0
    P_sat.pb = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.pb .* (sigma.apb + sigma.epb));
    
    % скорость вынужденного излучательного перехода для сигнала и ASE
    W21   = W21 + sum(psi.nase .* repmat(P(N.s+N.ase+N.pf+1: N.s+2*N.ase+N.pf, 1),1,size(psi.nase,2)) .*...
        repmat(sigma.ease',1,size(psi.nase,2)) ./ ((repmat(sigma.aase',1,size(psi.nase,2)) +...
        repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))),1);
    
    R24  = R24 + sum(psi.nase .* repmat(P(N.s+N.ase+N.pf+1: N.s+2*N.ase+N.pf, 1),1,size(psi.nase,2)) .*...
        repmat(sigma.esa_eta.*sigma.aase',1,size(psi.nase,2)) ./ (ph_const.h * ph_const.c * 10^25).* ...
        repmat(Lambda.ase',1,size(psi.nase,2)),1);

    for i = 1:length(Lambda.pb)
        if Lambda.pb(i) > 1450* 10^(-9)  % для накачки с длиной волны 1480 нм          
            
            % скорость вынужденного перехода с поглощением излучения для накачки
            R12  = R12 + psi.npb(i,:) .* P(N.s+2*N.ase+N.pf+i,1) .*  sigma.apb(i) ./ (sigma.apb(i) + sigma.epb(i)) ./ P_sat.pb(i);
            
            % скорость вынужденного излучательного перехода для накачки
            R21  = R21 + psi.npb(i,:) .* P(N.s+2*N.ase+N.pf+i,1) .*  sigma.epb(i) ./ (sigma.apb(i) + sigma.epb(i)) ./ P_sat.pb(i);
            
            % учет эффекта ESA для накачки на длине волны 1480 нм (для 980 нм должна быть большая мощность накачки)
            R24  = R24 + psi.npb(i,:) .* P(N.s+2*N.ase+N.pf+i,1) .* sigma.esa_eta*sigma.apb(i) ./ (ph_const.h * ph_const.c * 10^25) .* Lambda.pb(i);
        else                             % для накачки с длиной волны 980 нм
            
            % скорость вынужденного перехода с поглощением излучения для накачки
            R13   = R13 + psi.npb(i,:) .* P(N.s+2*N.ase+N.pf+i) .*  sigma.apb(i) ./ (sigma.apb(i) + sigma.epb(i)) ./ P_sat.pb(i);
            
            % скорость вынужденного излучательного перехода для накачки
            R31   = R31 + psi.npb(i,:) .* P(N.s+2*N.ase+N.pf+i) .*  sigma.epb(i) ./ (sigma.apb(i) + sigma.epb(i)) ./ P_sat.pb(i);
        end
    end
end

if isempty(Lambda.pb) == 1
    Lambda.pb = Lambda.pf;
end
%% Расчет населенностей четырехуровневой системы (расширение двухуровневой системы путем учета апконверсии и ESA)

% учет кооперативной апконверсии
C24 = (2.65 * n_sum * 0.1 + 3.38) * 10; % учет кооперативной апконверсии с уровня 2 на 4 (1-ое слагаемое:
                                        % коэффициент двухчастичной апконверсии, зависит от концентрации;
                                        % 2-ое слагаемое: коэффициент ветвления между переходом I(11/2)-I(15/2) (980 нм)
                                        % и переходом I(11/2)-I(13/2) (безызлучательный переход))

% для накачки с длиной волны 1480 нм 
if min(Lambda.pf) > 1450* 10^(-9) || min(Lambda.pb) > 1450* 10^(-9)
    
    % дискриминант для населенности n2
    D = (1 + W21./(W12+R12) + 1./(ph_const.tau2.*(W12+R12)) + R21./(W12+R12) + (ph_const.tau3+ph_const.tau4)*R24).^2 +...
    4.*(C24./(W12+R12) + (ph_const.tau3+ph_const.tau4)*C24);

    % массив населенностей уровня 2 : I(13/2), при разном положении в сердцевине
    n.second = ((-1).*(1 + W21./(W12+R12) + 1./(ph_const.tau2.*(W12+R12)) + R21./(W12+R12) + (ph_const.tau3+ph_const.tau4)*R24) + sqrt(D))...
    ./ (2*(C24./(W12+R12) + (ph_const.tau3+ph_const.tau4)*C24));

    % массив населенностей уровня 1 : I(15/2), при разном положении в сердцевине
    n.first = (n.second.*W21 + n.second/ph_const.tau2 + C24.*n.second.^2 + R21.*n.second) ./ (W12+R12);

    % массив населенностей уровня 4 : I(19/2), при разном положении в сердцевине
    n.fourth = ph_const.tau4 * (C24*n.second.^2 + n.second.*R24);

    % массив населенностей уровня 3 : I(11/2), при разном положении в сердцевине
    n.third = ph_const.tau3 / ph_const.tau4 .* n.fourth;

% для накачки с длиной волны 980 нм    
elseif max(Lambda.pf) < 1100* 10^(-9) || max(Lambda.pb) < 1100* 10^(-9) 
    
    % дискриминант для населенности n2
    D = (1/ph_const.tau2./(R13+W12) + W21./(R13+W12) + (R31./(R13+W12) + 1) * ph_const.tau3/ph_const.tau2 .* (R13 + W21*ph_const.tau2.*R13)...
        ./ (R13+W12+R31.*W12*ph_const.tau3) + 1).^2 + 4 * (C24./(R13+W12) + (R31./(R13+W12) + 1) * ph_const.tau3 * C24.*(2.*R13+W12)...
        ./ (R13+W12+R31.*W12*ph_const.tau3) + ph_const.tau4 * C24);
    
    % массив населенностей уровня 2 : I(13/2), при разном положении в сердцевине
    n.second = ((-1)*(1/ph_const.tau2./(R13+W12) + W21./(R13+W12) + (R31./(R13+W12) + 1)*ph_const.tau3/ph_const.tau2 .* (R13+W21*ph_const.tau2.*R13)...
        ./ (R13+W12+R31.*W12*ph_const.tau3) + 1) + sqrt(D))./2./(C24./(R13+W12) + (R31./(R13+W12) + 1)*ph_const.tau3 * C24.*(2.*R13+W12)...
        ./ (R13+W12+R31.*W12*ph_const.tau3) + ph_const.tau4 * C24);
    
    % массив населенностей уровня 4 : I(19/2), при разном положении в сердцевине
    n.fourth = C24 .* n.second.^2 * ph_const.tau4;
    
    % массив населенностей уровня 3 : I(11/2), при разном положении в сердцевине
    n.third  = ph_const.tau3/ph_const.tau2 .* (C24 .* n.second*ph_const.tau2.*(2.*R13+W12) + n.second.*R13 + W21.*n.second*ph_const.tau2.*R13)...
        ./ (R13+W12+R31.*W12*ph_const.tau3);
    
    % массив населенностей уровня 1 : I(15/2), при разном положении в сердцевине
    n.first = (n.third/ph_const.tau3 - n.fourth/ph_const.tau4 + R31.*n.third) ./ R13;

else
    % дискриминант для населенности n2
    D = (-R24./R13 + (1./R13/ph_const.tau3+R31./R13+1) * ph_const.tau3 .* ((R13+W12+R12).*R24 + R13/ph_const.tau2 + R13.*W21 + R13.*R21) ./...
        (R13+W12+R12+R31.*W12*ph_const.tau3+R31.*R12*ph_const.tau3) + 1 + ph_const.tau4.*R24).^2 + 4.*(-C24./R13+(1./R13/ph_const.tau3+R31./R13+1) *...
        ph_const.tau3 .* ((R13+W12+R12)*C24+R13*C24) ./ (R13+W12+R12+R31.*W12*ph_const.tau3+R31.*R12*ph_const.tau3) + ph_const.tau4*C24);
    
    % массив населенностей уровня 2 : I(13/2), при разном положении в сердцевине
    n.second = (-(-R24./R13+(1./R13/ph_const.tau3+R31./R13+1) * ph_const.tau3 .* ((R13+W12+R12).*R24+R13/ph_const.tau2+R13.*W21+R13.*R21) ./...
        (R13+W12+R12+R31.*W12*ph_const.tau3+R31.*R12*ph_const.tau3) + 1 + ph_const.tau4.*R24) + sqrt(D)) / 2. /(-C24./R13+(1./R13/ph_const.tau3+R31./R13+1) *...
        ph_const.tau3.*((R13+W12+R12)*C24+R13*C24) ./ (R13+W12+R12+R31.*W12*ph_const.tau3+R31.*R12*ph_const.tau3) + ph_const.tau4*C24);
    
    % массив населенностей уровня 4 : I(19/2), при разном положении в сердцевине
    n.fourth = ph_const.tau4 .* (C24 .* n.second.^2 + R24.*n.second);
    
    % массив населенностей уровня 3 : I(11/2), при разном положении в сердцевине
    n.third  = ph_const.tau3.*((C24.*n.second.^2+R24.*n.second).*(R13+W12+R12)+R13.*(n.second/ph_const.tau2+C24.*n.second.^2+W21.*n.second+R21.*n.second))...
        ./(R13+W12+R12+R31.*W12*ph_const.tau3+R31.*R12*ph_const.tau3);
    
    % массив населенностей уровня 1 : I(15/2), при разном положении в сердцевине
    n.first = (R21.*n.second + n.second/ph_const.tau2 - n.third/ph_const.tau3 + 2*C24.*n.second.^2 + W21.*n.second + R24 .*n.second)./(R12+W12);
end

% n2           = (W12 + R12) ./ (W12 + R12 + W21 + 1 / ph_const.tau + R21);                           % населенность верхнего уровня
% n1           = (1 / ph_const.tau + W21 + R21) ./ (W12 + R12 + W21 + 1 / ph_const.tau + R21);        % населенность нижнего уровня
end