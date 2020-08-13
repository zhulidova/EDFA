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

param         = sigma_lum_param();                                                      % ������ ������� ����������������� ����������
CrossLum      = 1.3 * param(:,3) .* exp((-1) * 10^(-21) * param(:,2) / ph_const.k / T); % ������ ������� �������������

% �����������
wl_smooth     = linspace(min(param(:,1)), max(param(:,1)), 80);                        % ����� ����� (80 ����� ���� ����)
avgCrossLum   = spline(param(:,1), CrossLum, wl_smooth);                               % ���������� ������ ������� �������������         

% ������ ������� ������� ����������
CrossAbs      = CrossLum .* exp((-1) * (ph_const.eV * 0.811 - ph_const.h * ph_const.c * 10^9 ./ param(:,1)) / (ph_const.k * T));
wl_smooth2    = linspace(min(param(:,1)), max(param(:,1)), 30);
avgCrossAbs   = spline(param(:,1), CrossAbs, wl_smooth2);                               % ���������� ������ �������       

% plot(wl_smooth,avgCrossLum);
% hold on;
% plot(wl_smooth2,avgCrossAbs);

sigma.ES     = interp1(wl_smooth, avgCrossLum, wl.S);       % ������� ������������� (������)  
sigma.AS     = interp1(wl_smooth2, avgCrossAbs, wl.S);      % ������� ���������� (������)
sigma.EPF    = interp1(wl_smooth, avgCrossLum, wl.PF);      % ������� ������������� (�������� �������)
sigma.APF    = interp1(wl_smooth2, avgCrossAbs, wl.PF);     % ������� ���������� (�������� �������)
sigma.EASE   = interp1(wl_smooth, avgCrossLum, wl.ASE);     % ������� ������������� (ASE)
sigma.AASE   = interp1(wl_smooth2, avgCrossAbs, wl.ASE);    % ������� ���������� (ASE)
sigma.EPB    = interp1(wl_smooth, avgCrossLum, wl.PB);      % ������� ������������� (��������� �������)
sigma.APB    = interp1(wl_smooth2, avgCrossAbs, wl.PB);     % ������� ���������� (��������� �������)

end