clear all; clc;
% P_s_sum         - суммарна€ мощность сигнала, дЅм
% P_in.S          - массив входных мощностей сигнала (суммарна€ мощность / кол-во длин волн), дЅм
% P_in.PF         - массив мощностей попутной накачки в 0, дЅм
% P_in.PB         - массив мощностей встречного ASE (на заданном массиве длин волн) в L, дЅм
%пустой массив, если накачка только попуна€
% r_edf           - радиус сердцевины активного волокна, м
% NA              - числова€ апертура активного волокна    
% n               - концентраци€ ионов эрби€, ppm
% L               - длина активного волокна, м
% T_c             - температура среды, ∞—
% low_gain_regime - механический выбор режима слабого сигнала (упрощенна€ система уравнений),
%0 - общий случай (без приближени€), 1 - режим слабого сигнала

% ‘ункции
% edfa_main       - основна€ расчетна€ функци€ EDFA
% smooth_res      - сглаживание результатов
% create_graf     - построение графиков 
P_p(1,:) = [1.85 3.89];
P_p(2,:) = [2.85 4.89];
P_p(3,:) = [3.85 5.89];
P_p(4,:) = [5.85 7.89];
P_p(5,:) = [7.85 9.89];
P_p(6,:) = [9.85 11.89];
P_p(7,:) = [11.85 13.89];
P_p(8,:) = [13.85 15.89];

Len = [6,7,8,9,10,11,12,13,14];
Temp = [-60, -50, -40, -20, 0, 25, 40, 70];
P_s = [-8 -3 2 7 10];
for k = 1:size(Temp,2)
for j = 1:size(P_p,1)
for i = 1:size(P_s,2)
    
P_s_sum = P_s(i);  % массив суммарных входных мощностей сигнала  
% параметры сигнала 
Lambda.s              = [1550.1	1551.7	1554.1	1557.33	1559.0	1560.58	1562.2	1564.65]* 10^(-9); % длины волн сигнала, нм
Pin.s                 = watt(P_s_sum) / length(Lambda.s) * ones(1,length(Lambda.s));        % входные мощности сигнала, дЅм

% параметры ASE 
Lambda_range          = [1528 1565] * 10^(-9);                           % спектр ASE, нм

% параметры накачки
Lambda.pf             = [1473 1480] * 10^(-9);                           % длины волн попутной накачки, нм
Pin.pf                = watt(P_p(j,:));

Lambda.pb             = [];                                               % длины волн встречной накачки, нм
                                                                         % (пустой массив, если накачка только попутна€)
Pin.pb                = [];                                              % входные мощности встречной накачки, дЅм

L    = 9;
N0.r                  = 25;
N0.ase                = 50;
N0.z                  = 50;
type                  = 'Pnppk_a2';

T_c = Temp(k);

%% ќсновна€ расчетна€ функци€
% [z, P, Gain, OSNR, NF, SPase]   = edfa(Pin, Lambda, Lambda_range, L, N0, type, T_c);
% Pout = dbm(P.s);


% сглаживание результатов моделировани€
% [wl_smooth.Model, G_smooth.Model, nf_smooth.Model] = smooth_res(Lambda.s', Gain, NF); 

% чтение экспериментальных данных
EXP   = exp_res(P_s_sum, T_c, dbm(sum(Pin.pf)),L);
  
[wl_smooth.Exp, G_smooth.Exp, nf_smooth.Exp]  = smooth_res(EXP(1,:), EXP(2,:), EXP(3,:));


%% оформление графиков
% create_graf(Lambda.s*10^9, Gain, "model", 'Gain', P_s_sum, Pin, 1, i); 
% hold on;
% create_graf(EXP(1,:), EXP(2,:), "exp", 'Gain', P_s_sum, Pin, 1, i, T_c); 
% hold on;
% create_graf(Lambda.s*10^9, NF, "model", 'NF', P_s_sum, Pin, 2, i); 
% hold on;
create_graf(wl_smooth.Exp, nf_smooth.Exp, "exp", 'NF', P_s_sum, Pin, 2, i, T_c);
hold on;
end
% Filename1 = "F_G" + sprintf("%.0f", T_c)+"_"+sprintf("%.0f", dbm(sum(Pin.pf)))+".emf";
Filename2 = "F_NF" + sprintf("%.0f", T_c)+"_"+sprintf("%.0f", dbm(sum(Pin.pf)))+".emf";
% saveas(gcf,Filename1);
saveas(gcf,Filename2);
hold off;
end
end
% plot(Len,G);
% hold on;
% plot(Len,G_exp);
% semilogx(watt(P_p)*10^3,T);
% hold on;
% P_in_exp = [-9.96 -7.99 -5.03 -2.97 0.04 2.02 4.97 6.99];
% P_out_exp = [-11.4400000000000 -9.28000000000000 -6.10000000000000 -3.85000000000000 -0.630000000000000 1.44000000000000 4.46000000000000 6.56000000000000];
% T_exp = watt(P_out_exp+0.33)./watt(P_in_exp);
% semilogx(watt(P_in_exp)*10^3,T_exp);
