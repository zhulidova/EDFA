%% ������� ������� ������� ������������� � ����������(�)
% �������� �� ������� ������� �������� �� ������ ���� �������
function sigma = sigma_lum(T, wl, ph_const)
% ����������
% param        - ������ ����������������� �������� ������� (param(:,1) - ����� �����, param(:,2)
%- ����������������� ����������)
% CrossTemp    - ������ ������� ������������� (param(:,1), CrossTemp)
% avgCrossTemp - ���������� ������ ������� ������������� 
% cross_abs    - ������ ������� ����������, ������������ �� ������� ���������� (param(:,1),
%cross_abs)
% sigma        - ��������� ������� ������������� (sigma.ES, sigma.EPF, sigma.EASE, sigma.EPB)
%� ���������� (sigma.AS, sigma.APF, sigma.AASE, sigma.APB)
param        = sigma_lum_param();
CrossTemp    = 1.3 * param(:,3) .* exp((-1) * 10^(-21) *param(:,2) / ph_const.k / T);

% �����������
x_smooth     = linspace(min(param(:,1)), max(param(:,1)), 80);                  % ����� ����� (500 ����� ���� ����)
avgCrossTemp = spline(param(:,1), CrossTemp, x_smooth);                        % ���������� ������ 1             

cross_abs    = avgCrossTemp .* exp((-1) * (ph_const.eV * 0.811 - ph_const.h * ph_const.c * 10^9 ./ (x_smooth)) / (ph_const.k * T));
% 
% plot(x_smooth,avgCrossTemp);
% hold on;
% plot(x_smooth,cross_abs);

sigma.ES     = interp1(x_smooth, avgCrossTemp, wl.S);       % ������� ������������� (������)  
sigma.AS     = interp1(x_smooth, cross_abs, wl.S);          % ������� ���������� (������)
sigma.EPF    = interp1(x_smooth, avgCrossTemp, wl.PF);      % ������� ������������� (�������� �������)
sigma.APF    = interp1(x_smooth, cross_abs, wl.PF);         % ������� ���������� (�������� �������)
sigma.EASE   = interp1(x_smooth, avgCrossTemp, wl.ASE);     % ������� ������������� (ASE)
sigma.AASE   = interp1(x_smooth, cross_abs, wl.ASE);        % ������� ���������� (ASE)
sigma.EPB    = interp1(x_smooth, avgCrossTemp, wl.PB);      % ������� ������������� (��������� �������)
sigma.APB    = interp1(x_smooth, cross_abs, wl.PB);         % ������� ���������� (��������� �������)

end