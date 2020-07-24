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
CrossTemp    = 2.5 * param(:,3) .* exp((-1) * 10^(-21) * param(:,2) / ph_const.k / T);

% �����������
a            = 4;                                             % �������� �����������            
coeff        = ones(1, a) / a;
avgCrossTemp = filter(coeff, 1, CrossTemp);
cross_abs    = avgCrossTemp .* exp((-1) * 10^4 * (1.602 * 0.852 - 6.626 ...
    * 300 ./ param(:,1)) / (1.38 * T));

% plot(param(:,1),[CrossTemp avgCrossTemp]);
% hold on;
% plot(param(:,1),cross_abs);

sigma.ES     = interp1(param(:,1), avgCrossTemp, wl.S);       % ������� ������������� (������)  
sigma.AS     = interp1(param(:,1), cross_abs, wl.S);          % ������� ���������� (������)
sigma.EPF    = interp1(param(:,1), avgCrossTemp, wl.PF);      % ������� ������������� (�������� �������)
sigma.APF    = interp1(param(:,1), cross_abs, wl.PF);         % ������� ���������� (�������� �������)
sigma.EASE   = interp1(param(:,1), avgCrossTemp, wl.ASE);     % ������� ������������� (ASE)
sigma.AASE   = interp1(param(:,1), cross_abs, wl.ASE);        % ������� ���������� (ASE)
sigma.EPB    = interp1(param(:,1), avgCrossTemp, wl.PB);      % ������� ������������� (��������� �������)
sigma.APB    = interp1(param(:,1), cross_abs, wl.PB);         % ������� ���������� (��������� �������)

end