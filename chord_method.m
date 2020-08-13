%% решение краевой задачи с помощью метода хорд (секущих, метод Ньютона)
% с помощью метода хорд находятся значения мощностей встречной накачки и встречного шума в начале линии
%(сведение задачи с граничными условиями к задаче Коши)объяснение метода - отчеты по модели EDFA/Жулидова. Отчет от 27.06)
function P_out = chord_method(L, wl, sigma, psi, N,w_edf, n_sum, P_in, low_gain_regime, ph_const, r_edf)
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
epsilon = 1e-6;
        ksi1(1,:)           = zeros(1, N.ASE);            % начальные условия для задачи Коши (ASE,правая часть итерации)
        ksi2(1,:)           = ones(1, N.ASE) * 10^(-9);   % начальные условия для задачи Коши (ASE,левая часть итерации)
        lamda1(1,:)         = zeros(1, N.PB);             % начальные условия для задачи Коши (P_p_backward,правая часть итерации)
        lamda2(1,:)         = ones(1, N.PB) * 0.0005;       % начальные условия для задачи Коши (P_p_backward,левая часть итерации)
        lamda               = ones(1, N.PB) * 1;
%% цикл для поиска недостающего условия для задачи Коши
     while abs(sum(lamda1) - sum(lamda)) > epsilon
         
         if low_gain_regime == 1 
            % решение упрощенной системы уравнений (режим слабого сигнала) 
            [z, P_out1]     = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, wl, sigma, N,w_edf, n_sum, ph_const, P_in), [0: 0.1: L], [undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF) ksi1 lamda1]);  
            [z, P_out2]     = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, wl, sigma, N,w_edf, n_sum, ph_const, P_in), [0: 0.1: L], [undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF) ksi2 lamda2]); 
         else 
            % решение подробной системы уравнений (общий случай)
            [z, P_out1]     = ode45(@(z,P) odu2(z, P, wl, sigma, psi, N, n_sum, ph_const, P_in, r_edf), [0: 0.1: L],[undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF) ksi1 lamda1]);  
            [z, P_out2]     = ode45(@(z,P) odu2(z, P, wl, sigma, psi, N, n_sum, ph_const, P_in, r_edf), [0: 0.1: L],[undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF) ksi2 lamda2]); 
         end
        
        % расчет значений встречной накачки и встречного ASE в конце линии
        P_ASEB1     = P_out1(size(P_out1,1), N.S + N.ASE + N.PF + 1 : N.S + N.ASE + N.PF + N.ASE);                           % полученное значение ASE в конце волокна, 1
        P_ASEB2     = P_out2(size(P_out2,1), N.S + N.ASE + N.PF + 1 : N.S + N.ASE + N.PF + N.ASE);                           % полученное значение ASE в конце волокна, 2
        P_PB1       = P_out1(size(P_out1,1), N.S + 2 * N.ASE + N.PF + 1 : N.S + N.ASE + N.PF + N.ASE + N.PB);                % полученное значение P_p_backward в конце волокна, 1
        P_PB2       = P_out2(size(P_out2,1), N.S + 2 * N.ASE + N.PF + 1 : N.S + N.ASE + N.PF + N.ASE + N.PB);                % полученное значение P_p_backward в конце волокна, 2
        
        % приближенное недостающее условие для задачи Коши
        ksi         = ksi1 + (ksi2 - ksi1) .* (undbm(P_in.ASEB) - P_ASEB1) ./ (P_ASEB2 - P_ASEB1);
        lamda       = lamda1 + (lamda2 - lamda1) .* (undbm(P_in.PB) - P_PB1) ./ (P_PB2 - P_PB1);
        m           = lamda1;
        lamda1      = lamda;
        lamda       = m;
        ksi1        = ksi;    
     end
     
     % решение полученной задачи Коши
     if low_gain_regime == 1
        % решение упрощенной системы уравнений (режим слабого сигнала)
        [z, P_out]  = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, wl, sigma, N,w_edf, n_sum, ph_const, P_in), [0: 0.1: L],[undbm(P_in.S)...
            undbm(P_in.ASEF) undbm(P_in.PF) ksi1 lamda1]);
        
     else
         % решение подробной системы уравнений (общий случай)
        [z, P_out]  = ode45(@(z,P) odu2(z, P, wl, sigma, psi, N, n_sum, ph_const, P_in, r_edf), [0: 0.1: L], [undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF) ksi1 lamda1]);
     end
end