%% ������� ���������� �������� 
% ���������� 
% color_map - ����� ������ ��� ��������
% P_s_sum   - ��������� �������� �������
% type("model" ��� "exp") - ���������� ��� �����
% spec("Gain" ��� "NF")   - ���������� ��� ������ � ������ ��������� � ������� ������������ ���
function [] = create_graf(Lambda, x,sigma, type, spec, P_s_sum, Pin, num, count, T_c)

color_map(1,:) = [1 0 0];
color_map(2,:) = [0 1 0];
color_map(3,:) = [0 0 1];
color_map(4,:) = [1 0.5 0];
color_map(5,:) = [1 0 0.5];
color_map(6,:) = [0.5 1 0];
color_map(7,:) = [0.5 0 1];
color_map(8,:) = [0 0.5 1];
color_map(9,:) = [0.5 0.5 1];

set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 14, 'DefaultAxesFontName', 'Arial');

% ������ - �������, ����������� - �������� ������

formatLegendPowerS = "P_{s}= " + sprintf("%.0f", P_s_sum) + " ���";
formatLegendPump = "";
if isempty(Pin.pb) == 0
    formatLegendPowerPB = sprintf("%.0f", round(dbm(sum(Pin.pb))));
    formatLegendPump = ", P_{p}^{b}= " + formatLegendPowerPB + " ��� ";
end
formatLegendPowerPF = sprintf("%.0f", round(dbm(sum(Pin.pf))));
formatLegendPump = formatLegendPump + " P_{p}= " + formatLegendPowerPF + " ��� ";
formatLegend1     = "P_{s}= ";
formatLegend2     = " ���";
formatLegendType1 = ", ������";
formatLegendType2 = ", �����������";
formatLegendFin1  = formatLegend1 + formatLegendPowerS + formatLegend2;
formatLegendFin2  = formatLegend1 + formatLegendPowerS + formatLegend2;

f = figure(num);
if type == "model"
    style                  = '-';
    if count > 1
        plot(Lambda, x, 'Color',color_map(count,:), 'LineStyle', style, 'LineWidth', 2);
    else
        plot(Lambda, x, 'Color',color_map(count,:), 'LineStyle', style, 'LineWidth', 2);
    end
elseif type == "exp"
    style                  = '-';
    e = errorbar(Lambda,x,sigma,'LineStyle','none','Color',color_map(count,:),'Marker','s','MarkerSize',10,'MarkerFaceColor',color_map(count,:),'CapSize',10);
end

hold on;
xlabel('����� �����, ��');
xlim([1529 1561]);
ylabel(spec + ", ��");
grid on;
% ����������� ������� � ��������� ��������� ������� � ��������� �������
main_title = "BROPA " + spec + "(" + formatLegendPowerS + "," + "T = " + sprintf("%.0f", T_c) + " �C" + ", EDF1)";
title(main_title); 
f.Position = [400 400 800 600];
set(gca,'position',[0.1 0.13 0.68 0.79])% �������� �������
end