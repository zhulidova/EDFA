%% ‘ункци€ построени€ графиков 
% ѕеременные 
% color_map - набор цветов дл€ спектров
% P_s_sum   - суммарна€ мощность сигнала
% type("model" или "exp") - определ€ет тип линии
% spec("Gain" или "NF")   - определ€ет тип кривой и делает назавание и подпись вертикальной оси
 function [] = create_graf(wl, x, type, spec, P_s_sum, P_in, num, count)

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

% составление легенды с суммарной мощностью сигнала и мощностью накачки
formatLegendPowerS = sprintf("%.0f", P_s_sum);
formatLegendPump = "";
if isempty(P_in.PB) == 0
    formatLegendPowerPB = sprintf("%.0f", round(dbm(sum(undbm(P_in.PB)))));
    formatLegendPump = ", P_{н}^{встр}= " + formatLegendPowerPB + " дЅм ";
end
    formatLegendPowerPF = sprintf("%.0f", round(dbm(sum(undbm(P_in.PF)))));
    formatLegendPump = formatLegendPump + ", P_{н}^{попут}= " + formatLegendPowerPF + " дЅм ";
formatLegend1     = "P_{с}^{вх}= ";
formatLegend2     = " дЅм";
formatLegendType1 = ", модель";
formatLegendType2 = ", эксперимент";
formatLegendFin1  = formatLegend1 + formatLegendPowerS + formatLegend2 + formatLegendPump + formatLegendType1 ;
formatLegendFin2  = formatLegend1 + formatLegendPowerS + formatLegend2 + formatLegendPump + formatLegendType2;
legend(formatLegendFin1, formatLegendFin2);
title(spec);                                                                        % название графика
end