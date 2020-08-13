%% ‘ункци€ построени€ графиков 
% ѕеременные 
% color_map - набор цветов дл€ спектров
% P_s_sum   - суммарна€ мощность сигнала
% type("model" или "exp") - определ€ет тип линии
% spec("Gain" или "NF")   - определ€ет тип кривой и делает назавание и подпись вертикальной оси
function [] = create_graf(wl, x, type, spec, P_s_sum, num, count)

color_map              = ['r', 'b', 'm', 'c'];
set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 14, 'DefaultTextFontName', 'Times New Roman');

% модель - пунктир, эксперимент - сплошна€ крива€
if type == "model"
style                  = '--';
else
style                  = '-';
end

figure(num);

plot(wl, x, color_map(count), 'LineStyle', style, 'LineWidth', 2);                        
hold on;                   
xlabel('ƒлина волны, нм'); 
xlim([1529 1561]);
ylabel(spec + ", дЅ");
grid on;

% составление легенды с суммарной мощностью сигнала
formatLegendPower = sprintf("%.0f",P_s_sum);             
formatLegend1     = "P_{с}^{вх}= ";
formatLegend2     = " дЅм";
formatLegendType1 = ", модель";
formatLegendType2 = ", эксперимент";
formatLegendFin1  = formatLegend1 + formatLegendPower + formatLegend2 + formatLegendType1;
formatLegendFin2  = formatLegend1 + formatLegendPower + formatLegend2 + formatLegendType2;
legend(formatLegendFin1, formatLegendFin2);
title(spec);                                                                        % название графика
end