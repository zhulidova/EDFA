function [G, nf, P_out] = edfa_main(P_in, wl, n, low_gain_regime, L, T_c, r_edf, NA, splices) 
%% основная расчетная функция
% Переменные
% ph_const      - структура с физическими константами
% T             - температура в кельвинах
% splices.wdm_p - потери на wdm для накачки (datasheet), дБм                         
% splices.wdm_s - потери на wdm для сигнала (datasheet), дБм                               
% splices.fiber - потери на сварке 
% n_sum         - концентрация в м^(-3)
% N             - структура с размерами массивов сигнала (N.S), накачек (N.PF и N.PB) и ASE (N.ASE)
% sigma         - структура с сечениями излучения сигнала (sigma.ES), накачки (sigma.EPF и sigma.EPB) и ASE (sigma.EASE)
%и сечениями поглощения сигнала (sigma.AS), накачки (sigma.APF и sigma.APB) и ASE (sigma.AASE)
% psi           - структура с функциями распередения интенсивности по радиусу волокна для сигнала (psi.S),
%накачки (psi.PF и psi.PB), ASE (psi.ASE) и нормированными функциями распередения (psi.NS, psi.PF, psi.PB, psi.ASE) 
% w_edf         - структура модовых радиусов для сигнала (w_edf.S), накачки (w_edf.PF и w_edf.PB) и ASE (w_edf.ASE) 

% Функции
% concentration(n)              - функция пересчета концентрации из ppm в м^(-3)
% sigma_lum(T, wl, N, ph_const) - функция расчета сечений поглощения и люминесценции для каждой длины волны 
%при определенной температуре
% bessel(r_edf, wl, N, NA)      - функция расчета функций распередения интенсивности по радиусу волокна для
%каждой дины волны (зависит от радиуса сердцевины и числовой апертуры активного волокна)
% splice_loss(P_in.S,..., NA)   - функция расчета реальных входных мощностей сигнала и накачки с учетом потерь
%на сварках и WDM
% odu2_unsaturated_gain_regime  - функция расчета значений производных для скоростных уравнений в приближении
%слабого сигнала 
% odu2                          - функция расчета значений производных для скоростных уравнений в общем случае 
% chord_method(L,..., ph_const) - функция с линейной интерполяцией начальных условий для сведения комбинированного
%случая или случая со встречной накачкой к задаче Коши 
% G(P_out, P_in.S)              - функция расчета коэффициента усиления
% nf(P_out, ph_const)           - функция расчета шум-фактора

ph_const.c             = 299792458;                                             % скорость света    
ph_const.h             = 6.626E-34;                                             % постоянная Планка
ph_const.tau           = 0.0102;                                                % время жизни на верхнем уровне
ph_const.k             = 1.38 * 10^(-23);                                       % постоянная Больцмана

% расчет параметров для ОДУ
T                      = 273.15 + T_c;                                          % температура в кельвинах
n_sum                  = concentration(n);                                      % концентрация в м^(-3)
N.S                    = length(wl.S);                                          % размер массива длин волн сигнала
N.PF                   = length(wl.PF);                                         % размер массива длин волн попутной накачки
N.ASE                  = length(wl.ASE);                                        % размер массива длин волн ASE
N.PB                   = length(wl.PB);                                         % размер массива длин встречной накачки
sigma                  = sigma_lum(T, wl, ph_const);                            % расчет сечений
[psi, w_edf]           = bessel(r_edf, wl, N, NA);                              % волновые функции интенсивности главной моды

    %% учет доп. потерь  

[P_in.S, P_in.PF] = splice_loss(P_in.S, P_in.PF, wl, r_edf, splices, NA);      % реальные входные мощности сигнала и накачки
                                                                              %с учетом потерь на сварках и WDM

%% решение сиситемы ОДУ
    tStart = cputime;                                                           % старт таймера
    if isempty(P_in.PB) == 1 
% случай только попуной накачки (задача Коши)
        if low_gain_regime == 1
% решение упрощенной системы уравнений (режим слабого сигнала)
            [z, P_out]     = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, wl, sigma, N, w_edf, n_sum, ph_const, P_in), [0: 0.1: L],[undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF)]);
        else
% решение подробной системы уравнений (общий случай)
            [z, P_out]     = ode45(@(z,P) odu2(z, P, wl, sigma, psi, N, n_sum, ph_const, P_in), [0: 0.1: L],[undbm(P_in.S) undbm(P_in.ASEF) undbm(P_in.PF)]);
        end
        G                  = Gain(P_out(size(P_out,1), 1 : N.S)', P_in.S);      % расчет КУ
        nf                 = NF(P_out(size(P_out,1), N.S + 1 : N.S + N.ASE),... % расчет ШФ
            G, wl, ph_const);                                                    
    else
% случай встречной накачки или комбинированный случай
        P_out              = chord_method(L, wl, sigma, psi, N, w_edf, n_sum, P_in, low_gain_regime, ph_const);
        G                  = Gain(P_out(size(P_out,1), 1 : N.S)', P_in.S);      % расчет КУ
        nf                 = NF(P_out(size(P_out,1), N.S + 1 : N.S + N.ASE),... % расчет ШФ
            G, wl, ph_const);                                                   
    end
    P_out = P_out / undb(splices.fiber) / undb(splices.wdm_s);      % реальные входные мощности сигнала и накачки
    time                   = cputime - tStart;                                  % конец таймера
    formatSpec             = 'Система ОДУ решена за %f секунд\n';                           
    fprintf(formatSpec,time)
end