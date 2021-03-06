%% ������� ���������� �������� 
% ���������� 
% color_map - ����� ������ ��� ��������
% P_s_sum   - ��������� �������� �������
% type("model" ��� "exp") - ���������� ��� �����
% spec("Gain" ��� "NF")   - ���������� ��� ������ � ������ ��������� � ������� ������������ ���
 function [] = create_graf(Lambda, x, type, spec, P_s_sum, Pin, num, count)

color_map              = ['r', 'b', 'm', 'c'];
set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 14, 'DefaultTextFontName', 'Times New Roman');

% ������ - �������, ����������� - �������� ������
if type == "model"
style                  = '--';
else
style                  = '-';
end

figure(num);

plot(Lambda, x, color_map(count), 'LineStyle', style, 'LineWidth', 2);                        
hold on;                   
xlabel('����� �����, ��'); 
xlim([1529 1561]);
ylabel(spec + ", ��");
grid on;

% ����������� ������� � ��������� ��������� ������� � ��������� �������
formatLegendPowerS = sprintf("%.0f", P_s_sum);
formatLegendPump = "";
if isempty(Pin.pb) == 0
    formatLegendPowerPB = sprintf("%.0f", round(dbm(sum(Pin.pb))));
    formatLegendPump = ", P_{�}^{����}= " + formatLegendPowerPB + " ��� ";
end
    formatLegendPowerPF = sprintf("%.0f", round(dbm(sum(Pin.pf))));
    formatLegendPump = formatLegendPump + ", P_{�}^{�����}= " + formatLegendPowerPF + " ��� ";
formatLegend1     = "P_{�}^{��}= ";
formatLegend2     = " ���";
formatLegendType1 = ", ������";
formatLegendType2 = ", �����������";
formatLegendFin1  = formatLegend1 + formatLegendPowerS + formatLegend2 + formatLegendPump + formatLegendType1 ;
formatLegendFin2  = formatLegend1 + formatLegendPowerS + formatLegend2 + formatLegendPump + formatLegendType2;
legend(formatLegendFin1, formatLegendFin2);
title(spec);                                                                        % �������� �������
end