%% ������� ����������� �������� (��� ��������)
function [x_smooth, y1_smooth, y2_smooth] = smooth_res(x,y1,y2)
x_smooth      = linspace(min(x), max(x), 500);                  % ����� ����� (500 ����� ���� ����)
y1_smooth     = spline(x, y1, x_smooth);                        % ���������� ������ 1              
y2_smooth     = spline(x, y2, x_smooth);                        % ���������� ������ 2
end