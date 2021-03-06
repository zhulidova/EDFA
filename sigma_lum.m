%% ������� ������� ������� ������������� � ����������(�)
% �������� �� ������� ������� �������� �� ������ ���� �������
function sigma = sigma_lum(T, Lambda, ph_const)
% ����������
% param        - ������ ����������������� �������� ������� (param(:,1) - ����� �����, param(:,2)
%- ����������������� ����������)
% CrossTemp    - ������ ������� ������������� (param(:,1), CrossTemp)
% avgCrossTemp - ���������� ������ ������� ������������� 
% cross_abs    - ������ ������� ����������, ������������ �� ������� ���������� (param(:,1),
%cross_abs)
% sigma        - ��������� ������� ������������� (sigma.ES, sigma.EPF, sigma.EASE, sigma.EPB)
%� ���������� (sigma.AS, sigma.APF, sigma.AASE, sigma.APB)

param         = sigma_lum_param();                                                      % ������ ������� ����������������� ����������
CrossLum      = 2.35 * param(:,3) .* exp((-1) * 10^(-21) * param(:,2) / ph_const.k / T); % ������ ������� �������������

% �����������
wl_smooth     = linspace(min(param(:,1)), max(param(:,1)), 100);                        % ����� ����� (80 ����� ���� ����)
avgCrossLum   = spline(param(:,1), CrossLum, wl_smooth);                               % ���������� ������ ������� �������������         

% ������ ������� ������� ����������
CrossAbs      = avgCrossLum .* exp((-1) * (ph_const.eV * 0.809 - ph_const.h * ph_const.c * 10^9 ./ wl_smooth) / (ph_const.k * T));
wl_smooth2    = linspace(min(param(:,1)), max(param(:,1)), 100);
avgCrossAbs   = spline(wl_smooth, CrossAbs, wl_smooth2);                               % ���������� ������ �������       

% plot(wl_smooth,avgCrossLum);
% hold on;
% plot(wl_smooth2,avgCrossAbs);

sigma.es     = interp1(wl_smooth, avgCrossLum, Lambda.s * 10^(9));       % ������� ������������� (������)  
sigma.as     = interp1(wl_smooth2, avgCrossAbs, Lambda.s * 10^(9));      % ������� ���������� (������)
sigma.epf    = interp1(wl_smooth, avgCrossLum, Lambda.pf * 10^(9));      % ������� ������������� (�������� �������)
sigma.apf    = interp1(wl_smooth2, avgCrossAbs, Lambda.pf * 10^(9));     % ������� ���������� (�������� �������)
sigma.ease   = interp1(wl_smooth, avgCrossLum, Lambda.ase * 10^(9));     % ������� ������������� (ASE)
sigma.aase   = interp1(wl_smooth2, avgCrossAbs, Lambda.ase * 10^(9));    % ������� ���������� (ASE)
sigma.epb    = interp1(wl_smooth, avgCrossLum, Lambda.pb * 10^(9));      % ������� ������������� (��������� �������)
sigma.apb    = interp1(wl_smooth2, avgCrossAbs, Lambda.pb * 10^(9));     % ������� ���������� (��������� �������)

end