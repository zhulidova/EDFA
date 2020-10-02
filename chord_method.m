%% решение краевой задачи с помощью метода хорд (секущих, метод Ньютона)
% с помощью метода хорд находятся значения мощностей встречной накачки и встречного шума в начале линии
%(сведение задачи с граничными условиями к задаче Коши)объяснение метода - отчеты по модели EDFA/Жулидова. Отчет от 27.06)
function [z, P_out] = chord_method(L, Lambda, sigma, psi, N,w_edf, n_sum, Pin, low_gain_regime, ph_const, r_edf, N0)
% Переменные 
% epsilon           - значение точности итерационного метода
% ksi1, ksi2        - значения начальных условий для ASE
% lamda1, lamda2    - значения начальных условий для накачки
% lamda             - предыдущая итерация lamda1 
% P_ASEB1 и P_ASEB2 - значения мощности встречных шумов в конце линии
% P_PB1 и P_PB2     - значения мощности встречных накачек в конце линии

% Функции
% odu2_unsaturated_gain_regime  - функция расчета значений производных для скоростных уравнений в приближении
%слабого сигнала 
% odu2                          - функция расчета значений производных для скоростных уравнений в общем случае
epsilon = 1e-9;
        ksi1(1,:)           = zeros(1, N.ase);            % начальные условия для задачи Коши (ASE,правая часть итерации)
        ksi2(1,:)           = ones(1, N.ase) * 10^(-8);   % начальные условия для задачи Коши (ASE,левая часть итерации)
        lamda1(1,:)         = zeros(1, N.pb)* 0.001;             % начальные условия для задачи Коши (P_p_backward,правая часть итерации)
        lamda2(1,:)         = ones(1, N.pb) * 0.005;     % начальные условия для задачи Коши (P_p_backward,левая часть итерации)
        lamda               = ones(1, N.pb) * 1;
%% цикл для поиска недостающего условия для задачи Коши
     while abs(sum(lamda1) - sum(lamda)) > epsilon
         if lamda1(1,1)< 0 || lamda1(1,2) < 0
             lamda1(1,:)         = ones(1, N.pb) * 3*10^(-6);             % начальные условия для задачи Коши (P_p_backward,правая часть итерации)
             lamda2(1,:)         = ones(1, N.pb) * 2*10^(-5);
         end
         if low_gain_regime == 1 
            % решение упрощенной системы уравнений (режим слабого сигнала) 
            [z, P_out1]     = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, Lambda, sigma, N,w_edf, n_sum, ph_const, Pin), [0: L/N0.z: L], [Pin.s Pin.asef Pin.pf ksi1 lamda1]);  
            [z, P_out2]     = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, Lambda, sigma, N,w_edf, n_sum, ph_const, Pin), [0: L/N0.z: L], [Pin.s Pin.asef Pin.pf ksi2 lamda2]); 
         else 
            % решение подробной системы уравнений (общий случай)
            [z, P_out1]     = ode45(@(z,P) odu2(z, P, Lambda, sigma, psi, N, n_sum, ph_const, r_edf, N0), [0: L/N0.z: L],[Pin.s Pin.asef Pin.pf ksi1 lamda1]);  
            [z, P_out2]     = ode45(@(z,P) odu2(z, P, Lambda, sigma, psi, N, n_sum, ph_const, r_edf, N0), [0: L/N0.z: L],[Pin.s Pin.asef Pin.pf ksi2 lamda2]); 
         end
        
        % расчет значений встречной накачки и встречного ASE в конце линии
        P_ASEB1     = P_out1(size(P_out1,1), N.s+N.ase+N.pf+1: N.s+N.ase+N.pf+N.ase);                       % полученное значение ASE в конце волокна, 1
        P_ASEB2     = P_out2(size(P_out2,1), N.s+N.ase+N.pf+1: N.s+N.ase+N.pf+N.ase);                       % полученное значение ASE в конце волокна, 2
        P_PB1       = P_out1(size(P_out1,1), N.s+2*N.ase+N.pf+1: N.s+N.ase+N.pf+N.ase+N.pb);                % полученное значение P_p_backward в конце волокна, 1
        P_PB2       = P_out2(size(P_out2,1), N.s+2*N.ase+N.pf+1: N.s+N.ase+N.pf+N.ase+N.pb);                % полученное значение P_p_backward в конце волокна, 2
        
        % приближенное недостающее условие для задачи Коши
        ksi         = ksi1 + (ksi2 - ksi1) .* (Pin.aseb - P_ASEB1) ./ (P_ASEB2 - P_ASEB1);
        lamda       = lamda1 + (lamda2 - lamda1) .* (Pin.pb - P_PB1) ./ (P_PB2 - P_PB1);
        m           = lamda1;
        lamda1      = lamda;
        lamda       = m;
        ksi1        = ksi;    
     end
     
     % решение полученной задачи Коши
     if low_gain_regime == 1
        % решение упрощенной системы уравнений (режим слабого сигнала)
        [z, P_out]  = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, Lambda, sigma, N,w_edf, n_sum, ph_const, Pin), [0: L/N0.z: L],[Pin.s Pin.asef Pin.pf ksi1 lamda1]);
        
     else
         % решение подробной системы уравнений (общий случай)
        [z, P_out]  = ode45(@(z,P) odu2(z, P, Lambda, sigma, psi, N, n_sum, ph_const, r_edf, N0), [0: L/N0.z: L], [Pin.s Pin.asef Pin.pf ksi1 lamda1]);
     end
end