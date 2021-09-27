%% ‘ункци€ построени€ графиков 
% ѕеременные 
% color_map - набор цветов дл€ спектров
% P_s_sum   - суммарна€ мощность сигнала
% type("model" или "exp") - определ€ет тип линии
% spec("Gain" или "NF")   - определ€ет тип кривой и делает назавание и подпись вертикальной оси
function [] = create_graf(Lambda, x, type, spec, P_s_sum, Pin, num, count, T_c)

color_map(1,:) = [1 0 0];
color_map(2,:) = [0 1 0];
color_map(3,:) = [0 0 1];
color_map(4,:) = [1 0.5 0];
color_map(5,:) = [1 0 0.5];
color_map(6,:) = [0.5 1 0];
color_map(7,:) = [0.5 0 1];
color_map(8,:) = [0 0.5 1];
color_map(9,:) = [0.5 0.5 1];

set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 14, 'DefaultTextFontName', 'Times New Roman');

% модель - пунктир, эксперимент - сплошна€ крива€

formatLegendPowerS = sprintf("%.0f", P_s_sum);
formatLegendPump = "";
if isempty(Pin.pb) == 0
    formatLegendPowerPB = sprintf("%.0f", round(dbm(sum(Pin.pb))));
    formatLegendPump = ", P_{p}^{b}= " + formatLegendPowerPB + " дЅм ";
end
formatLegendPowerPF = sprintf("%.0f", round(dbm(sum(Pin.pf))));
formatLegendPump = formatLegendPump + " P_{p}= " + formatLegendPowerPF + " дЅм ";
formatLegend1     = "P_{s}= ";
formatLegend2     = " дЅм";
formatLegendType1 = ", модель";
formatLegendType2 = ", эксперимент";
formatLegendFin1  = formatLegend1 + formatLegendPowerS + formatLegend2;
formatLegendFin2  = formatLegend1 + formatLegendPowerS + formatLegend2;

figure(num);
if type == "model"
    style                  = '--';
    if count > 1
        plot(Lambda, x, 'Color',color_map(count,:), 'LineStyle', style, 'LineWidth', 2);
    else
        plot(Lambda, x, 'Color',color_map(count,:), 'LineStyle', style, 'LineWidth', 2,'DisplayName',formatLegendFin1);
    end
elseif type == "exp"
    style                  = '-';
    plot(Lambda, x, 'Color',color_map(count,:), 'LineStyle', style, 'LineWidth', 2,'DisplayName',formatLegendFin2);
end
legend('Location', [0.8 0.71 0.15 0.05]);
hold on;
xlabel('ƒлина волны, нм');
xlim([1550 1565]);
ylabel(spec + ", дЅ");
grid on;
% составление легенды с суммарной мощностью сигнала и мощностью накачки
main_title = spec + " FROPA" + "(" + formatLegendPump + "," + "T = " + sprintf("%.0f", T_c) + " ∞C" + ")";
title(main_title); 

set(gca,'position',[0.1 0.13 0.65 0.78])% название графика
end