%% функция расчета сечения люминесценции и поглощения(Т)
% выбирает из спектра сечений значения на длинах волн сигнала
function sigma = sigma_lum(T, wl, ph_const)
% Переменные
% param        - массив экспериментальных значений сечений (param(:,1) - длина волны, param(:,2)
%- экспериментальные постоянные)
% CrossTemp    - спектр сечения люминесценции (param(:,1), CrossTemp)
% avgCrossTemp - сглаженный спектр сечения люминесценции 
% cross_abs    - спектр сечения поглощения, рассчитанный по формуле Маккамбера (param(:,1),
%cross_abs)
% sigma        - структура сечений люминесценции (sigma.ES, sigma.EPF, sigma.EASE, sigma.EPB)
%и поглощения (sigma.AS, sigma.APF, sigma.AASE, sigma.APB)
param        = sigma_lum_param();
CrossTemp    = 1.3 * param(:,3) .* exp((-1) * 10^(-21) *param(:,2) / ph_const.k / T);

% сглаживание
x_smooth     = linspace(min(param(:,1)), max(param(:,1)), 80);                  % новая сетка (500 точек длин волн)
avgCrossTemp = spline(param(:,1), CrossTemp, x_smooth);                        % сглаженный спектр 1             

cross_abs    = avgCrossTemp .* exp((-1) * (ph_const.eV * 0.811 - ph_const.h * ph_const.c * 10^9 ./ (x_smooth)) / (ph_const.k * T));
% 
% plot(x_smooth,avgCrossTemp);
% hold on;
% plot(x_smooth,cross_abs);

sigma.ES     = interp1(x_smooth, avgCrossTemp, wl.S);       % сечение люминесценции (сигнал)  
sigma.AS     = interp1(x_smooth, cross_abs, wl.S);          % сечение поглощения (сигнал)
sigma.EPF    = interp1(x_smooth, avgCrossTemp, wl.PF);      % сечение люминесценции (попутная накачка)
sigma.APF    = interp1(x_smooth, cross_abs, wl.PF);         % сечение поглощения (попутная накачка)
sigma.EASE   = interp1(x_smooth, avgCrossTemp, wl.ASE);     % сечение люминесценции (ASE)
sigma.AASE   = interp1(x_smooth, cross_abs, wl.ASE);        % сечение поглощения (ASE)
sigma.EPB    = interp1(x_smooth, avgCrossTemp, wl.PB);      % сечение люминесценции (встречная накачка)
sigma.APB    = interp1(x_smooth, cross_abs, wl.PB);         % сечение поглощения (встречная накачка)

end