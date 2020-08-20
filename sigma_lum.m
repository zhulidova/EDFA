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

param         = sigma_lum_param();                                                      % чтение массива экспериментальных параметров
CrossLum      = 2.5 * param(:,3) .* exp((-1) * 10^(-21) * param(:,2) / ph_const.k / T); % расчет сечения люминесценции

% сглаживание
wl_smooth     = linspace(min(param(:,1)), max(param(:,1)), 100);                        % новая сетка (80 точек длин волн)
avgCrossLum   = spline(param(:,1), CrossLum, wl_smooth);                               % сглаженный спектр сечения люминесценции         

% расчет спектра сечения поглощения
CrossAbs      = avgCrossLum .* exp((-1) * (ph_const.eV * 0.809 - ph_const.h * ph_const.c * 10^9 ./ wl_smooth) / (ph_const.k * T));
wl_smooth2    = linspace(min(param(:,1)), max(param(:,1)), 100);
avgCrossAbs   = spline(wl_smooth, CrossAbs, wl_smooth2);                               % сглаженный спектр сечения       

% plot(wl_smooth,avgCrossLum);
% hold on;
% plot(wl_smooth2,avgCrossAbs);

sigma.ES     = interp1(wl_smooth, avgCrossLum, wl.S);       % сечение люминесценции (сигнал)  
sigma.AS     = interp1(wl_smooth2, avgCrossAbs, wl.S);      % сечение поглощения (сигнал)
sigma.EPF    = interp1(wl_smooth, avgCrossLum, wl.PF);      % сечение люминесценции (попутная накачка)
sigma.APF    = interp1(wl_smooth2, avgCrossAbs, wl.PF);     % сечение поглощения (попутная накачка)
sigma.EASE   = interp1(wl_smooth, avgCrossLum, wl.ASE);     % сечение люминесценции (ASE)
sigma.AASE   = interp1(wl_smooth2, avgCrossAbs, wl.ASE);    % сечение поглощения (ASE)
sigma.EPB    = interp1(wl_smooth, avgCrossLum, wl.PB);      % сечение люминесценции (встречная накачка)
sigma.APB    = interp1(wl_smooth2, avgCrossAbs, wl.PB);     % сечение поглощения (встречная накачка)

end