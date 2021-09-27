function sigma = sigma_lum(T, Lambda, ph_const,n_sum)
%% Переменные
% sigma_em     - массив экспериментальных значений сечений (sigma_em(:,1) - длина волны, sigma_em(:,2),
% sigma_em(:,3), sigma_em(:,4), sigma_em(:,5), sigma_em(:,6),sigma_em(:,7),- экспериментальные сечения 
% люминесценции при температурах -40, -20, 0, 25, 40, 70 соответственно)
% sigma_abs     - массив экспериментальных значений сечений (sigma_abs(:,1) - длина волны, sigma_abs(:,2),
% sigma_abs(:,3), sigma_abs(:,4), sigma_abs(:,5), sigma_abs(:,6),sigma_abs(:,7),- экспериментальные сечения 
% поглощения при температурах -40, -20, 0, 25, 40, 70 соответственно)
% CrossLum     - спектр сечения люминесценции при заданной температуре
% CrossAbs     - спектр сечения поглощения при заданной температуре
% sigma        - структура сечений люминесценции (sigma.es, sigma.epf, sigma.ease, sigma.epb, sigma.esapf, sigma.esapb)
% и поглощения (sigma.as, sigma.apf, sigma.aase, sigma.apb)


sigma_em       = sigma_em_temp();  % чтение экспериментальных данных
sigma_abs      = sigma_abs_temp();

param_em(:,1)  = sigma_em(:,1);    % экспериментальная длина волны
param_abs(:,1) = sigma_abs(:,1);

% выбор двух близлежащих экспериментальных спектров и расчет параметров для сечения при T  
if T <= -20 + 273.15
    x = log(sigma_em(:,2)./sigma_em(:,3));
    param_em(:,2)  = 10^21*x * ph_const.k/(1/233.15 - 1/253.15);
    param_em(:,3)  = sigma_em(:,2) .* exp(-(1).*param_em(:,2)*10^(-21)/ph_const.k/233.15);
    y = log(sigma_abs(:,2)./sigma_abs(:,3));
    param_abs(:,2) = 10^21*y * ph_const.k/(1/233.15 - 1/253.15);
    param_abs(:,3) = sigma_abs(:,2) .* exp(-(1).*param_abs(:,2)*10^(-21)/ph_const.k/233.15);
elseif T > -20 + 273.15 && T <= 0 + 273.15
    x = log(sigma_em(:,3)./sigma_em(:,4));
    param_em(:,2)  = 10^21*x * ph_const.k/(1/253.15 - 1/273.15);
    param_em(:,3)  = sigma_em(:,3) .* exp(-(1).*param_em(:,2)*10^(-21)/ph_const.k/253.15);
    y = log(sigma_abs(:,3)./sigma_abs(:,4));
    param_abs(:,2) = 10^21*y * ph_const.k/(1/253.15 - 1/273.15);
    param_abs(:,3) = sigma_abs(:,3) .* exp(-(1).*param_abs(:,2)*10^(-21)/ph_const.k/253.15);
elseif T > 0 + 273.15 && T <= 25 + 273.15
    x = log(sigma_em(:,4)./sigma_em(:,5));
    param_em(:,2)  = 10^21*x * ph_const.k/(1/273.15 - 1/298.15);
    param_em(:,3)  = sigma_em(:,4) .* exp(-(1).*param_em(:,2)*10^(-21)/ph_const.k/273.15);
    y = log(sigma_abs(:,4)./sigma_abs(:,5));
    param_abs(:,2) = 10^21*y * ph_const.k/(1/273.15 - 1/298.15);
    param_abs(:,3) = sigma_abs(:,4) .* exp(-(1).*param_abs(:,2)*10^(-21)/ph_const.k/273.15);
elseif T > 25 + 273.15 && T <= 40 + 273.15
    x = log(sigma_em(:,5)./sigma_em(:,6));
    param_em(:,2)  = 10^21*x * ph_const.k/(1/298.15 - 1/313.15);
    param_em(:,3)  = sigma_em(:,5) .* exp(-(1).*param_em(:,2)*10^(-21)/ph_const.k/298.15);
    y = log(sigma_abs(:,5)./sigma_abs(:,6));
    param_abs(:,2) = 10^21*y * ph_const.k/(1/298.15 - 1/313.15);
    param_abs(:,3) = sigma_abs(:,5) .* exp(-(1).*param_abs(:,2)*10^(-21)/ph_const.k/298.15);
elseif T > 40 + 273.15
    x = log(sigma_em(:,6)./sigma_em(:,7));
    param_em(:,2)  = 10^21*x * ph_const.k/(1/313.15 - 1/343.15);
    param_em(:,3)  = sigma_em(:,6) .* exp(-(1).*param_em(:,2)*10^(-21)/ph_const.k/313.15);
    y = log(sigma_abs(:,6)./sigma_abs(:,7));
    param_abs(:,2) = 10^21*y * ph_const.k/(1/313.15 - 1/343.15);
    param_abs(:,3) = sigma_abs(:,6) .* exp(-(1).*param_abs(:,2)*10^(-21)/ph_const.k/313.15);
end

CrossLum       = 1/1.2*param_em(:,3) .* exp(10^(-21) * param_em(:,2) / ph_const.k / T); % расчет сечения люминесценции

CrossAbs       = 1/3.7*param_abs(:,3) .* exp(10^(-21) * param_abs(:,2) / ph_const.k / T); % расчет сечения люминесценции
% plot(param_abs(:,1),CrossAbs);
% hold on;
% plot(param_em(:,1),CrossLum);

sigma.es       = interp1(param_em(:,1), CrossLum, Lambda.s * 10^(9));             % сечение люминесценции (сигнал)
sigma.as       = interp1(param_abs(:,1), CrossAbs, Lambda.s * 10^(9));            % сечение поглощения (сигнал)
sigma.apf      = interp1(param_abs(:,1), CrossAbs, Lambda.pf * 10^(9));           % сечение поглощения (попутная накачка)
sigma.ease     = interp1(param_em(:,1), CrossLum, Lambda.ase * 10^(9));           % сечение люминесценции (ASE)
sigma.aase     = interp1(param_abs(:,1), CrossAbs, Lambda.ase * 10^(9));          % сечение поглощения (ASE)
sigma.apb      = interp1(param_abs(:,1), CrossAbs, Lambda.pb * 10^(9));           % сечение поглощения (встречная накачка)
sigma.esa_eta    = 0.1;                   % отношение сечения перехода с 2 на 4 уровень к сечению перехода с 1 на 2 уровень, эффект ESA
                                                                                  % для попутной накачки 
sigma.bg       = 6.31 * 10^(-3) / 10/log10(exp(1))/n_sum;                          

%% условие для накачек на 1480 и 980 нм
for i = 1:length(Lambda.pf)
    if Lambda.pf(i) > 1450* 10^(-9)
        sigma.epf(i)  = interp1(param_em(:,1), CrossLum, Lambda.pf(i) * 10^(9));            % сечение люминесценции (попутная накачка)
    else
        sigma.epf(i)  = 0;
    end
end
if isempty(Lambda.pb) == 0
    for i = 1:length(Lambda.pb)
        if min(Lambda.pb) > 1450* 10^(-9)
            sigma.epb(i)  = interp1(param_em(:,1), CrossLum, Lambda.pb(i) * 10^(9));            % сечение люминесценции (встречная накачка)
        else
            sigma.epb(i)  = 0;
        end
    end
end
    
