function sigma = sigma_lum(T, Lambda, ph_const,n_sum)
%% ����������
% sigma_em     - ������ ����������������� �������� ������� (sigma_em(:,1) - ����� �����, sigma_em(:,2),
% sigma_em(:,3), sigma_em(:,4), sigma_em(:,5), sigma_em(:,6),sigma_em(:,7),- ����������������� ������� 
% ������������� ��� ������������ -40, -20, 0, 25, 40, 70 ��������������)
% sigma_abs     - ������ ����������������� �������� ������� (sigma_abs(:,1) - ����� �����, sigma_abs(:,2),
% sigma_abs(:,3), sigma_abs(:,4), sigma_abs(:,5), sigma_abs(:,6),sigma_abs(:,7),- ����������������� ������� 
% ���������� ��� ������������ -40, -20, 0, 25, 40, 70 ��������������)
% CrossLum     - ������ ������� ������������� ��� �������� �����������
% CrossAbs     - ������ ������� ���������� ��� �������� �����������
% sigma        - ��������� ������� ������������� (sigma.es, sigma.epf, sigma.ease, sigma.epb, sigma.esapf, sigma.esapb)
% � ���������� (sigma.as, sigma.apf, sigma.aase, sigma.apb)


sigma_em       = sigma_em_temp();  % ������ ����������������� ������
sigma_abs      = sigma_abs_temp();

param_em(:,1)  = sigma_em(:,1);    % ����������������� ����� �����
param_abs(:,1) = sigma_abs(:,1);

% ����� ���� ����������� ����������������� �������� � ������ ���������� ��� ������� ��� T  
if T <= -20 + 273.15
    x = log(sigma_em(:,2)./sigma_em(:,3));
    param_em(:,2)  = 10^21*x * ph_const.k/(1/233.15 - 1/253.15);
    param_em(:,3)  = sigma_em(:,2) .* exp(-(1).*param_em(:,2)*10^(-21)/ph_const.k/233.15);
    y = log(sigma_abs(:,2)./sigma_abs(:,3));
    param_abs(:,2) = 10^21*y * ph_const.k/(1/233.15 - 1/253.15);
    param_abs(:,3) = sigma_abs(:,2) .* exp(-(1).*param_abs(:,2)*10^(-21)/ph_const.k/233.15);
elseif T > -20 + 273.15 && T <= 0 + 273.15
    x = log(sigma_em(:,3)./sigma_em(:,4));
    param_em(:,2)  = 10^21*x * ph_const.k/(1/253.15 - 1/273.15);
    param_em(:,3)  = sigma_em(:,3) .* exp(-(1).*param_em(:,2)*10^(-21)/ph_const.k/253.15);
    y = log(sigma_abs(:,3)./sigma_abs(:,4));
    param_abs(:,2) = 10^21*y * ph_const.k/(1/253.15 - 1/273.15);
    param_abs(:,3) = sigma_abs(:,3) .* exp(-(1).*param_abs(:,2)*10^(-21)/ph_const.k/253.15);
elseif T > 0 + 273.15 && T <= 25 + 273.15
    x = log(sigma_em(:,4)./sigma_em(:,5));
    param_em(:,2)  = 10^21*x * ph_const.k/(1/273.15 - 1/298.15);
    param_em(:,3)  = sigma_em(:,4) .* exp(-(1).*param_em(:,2)*10^(-21)/ph_const.k/273.15);
    y = log(sigma_abs(:,4)./sigma_abs(:,5));
    param_abs(:,2) = 10^21*y * ph_const.k/(1/273.15 - 1/298.15);
    param_abs(:,3) = sigma_abs(:,4) .* exp(-(1).*param_abs(:,2)*10^(-21)/ph_const.k/273.15);
elseif T > 25 + 273.15 && T <= 40 + 273.15
    x = log(sigma_em(:,5)./sigma_em(:,6));
    param_em(:,2)  = 10^21*x * ph_const.k/(1/298.15 - 1/313.15);
    param_em(:,3)  = sigma_em(:,5) .* exp(-(1).*param_em(:,2)*10^(-21)/ph_const.k/298.15);
    y = log(sigma_abs(:,5)./sigma_abs(:,6));
    param_abs(:,2) = 10^21*y * ph_const.k/(1/298.15 - 1/313.15);
    param_abs(:,3) = sigma_abs(:,5) .* exp(-(1).*param_abs(:,2)*10^(-21)/ph_const.k/298.15);
elseif T > 40 + 273.15
    x = log(sigma_em(:,6)./sigma_em(:,7));
    param_em(:,2)  = 10^21*x * ph_const.k/(1/313.15 - 1/343.15);
    param_em(:,3)  = sigma_em(:,6) .* exp(-(1).*param_em(:,2)*10^(-21)/ph_const.k/313.15);
    y = log(sigma_abs(:,6)./sigma_abs(:,7));
    param_abs(:,2) = 10^21*y * ph_const.k/(1/313.15 - 1/343.15);
    param_abs(:,3) = sigma_abs(:,6) .* exp(-(1).*param_abs(:,2)*10^(-21)/ph_const.k/313.15);
end

CrossLum       = 1/1.2*param_em(:,3) .* exp(10^(-21) * param_em(:,2) / ph_const.k / T); % ������ ������� �������������

CrossAbs       = 1/3.7*param_abs(:,3) .* exp(10^(-21) * param_abs(:,2) / ph_const.k / T); % ������ ������� �������������
% plot(param_abs(:,1),CrossAbs);
% hold on;
% plot(param_em(:,1),CrossLum);

sigma.es       = interp1(param_em(:,1), CrossLum, Lambda.s * 10^(9));             % ������� ������������� (������)
sigma.as       = interp1(param_abs(:,1), CrossAbs, Lambda.s * 10^(9));            % ������� ���������� (������)
sigma.apf      = interp1(param_abs(:,1), CrossAbs, Lambda.pf * 10^(9));           % ������� ���������� (�������� �������)
sigma.ease     = interp1(param_em(:,1), CrossLum, Lambda.ase * 10^(9));           % ������� ������������� (ASE)
sigma.aase     = interp1(param_abs(:,1), CrossAbs, Lambda.ase * 10^(9));          % ������� ���������� (ASE)
sigma.apb      = interp1(param_abs(:,1), CrossAbs, Lambda.pb * 10^(9));           % ������� ���������� (��������� �������)
sigma.esa_eta    = 0.1;                   % ��������� ������� �������� � 2 �� 4 ������� � ������� �������� � 1 �� 2 �������, ������ ESA
                                                                                  % ��� �������� ������� 
sigma.bg       = 6.31 * 10^(-3) / 10/log10(exp(1))/n_sum;                          

%% ������� ��� ������� �� 1480 � 980 ��
for i = 1:length(Lambda.pf)
    if Lambda.pf(i) > 1450* 10^(-9)
        sigma.epf(i)  = interp1(param_em(:,1), CrossLum, Lambda.pf(i) * 10^(9));            % ������� ������������� (�������� �������)
    else
        sigma.epf(i)  = 0;
    end
end
if isempty(Lambda.pb) == 0
    for i = 1:length(Lambda.pb)
        if min(Lambda.pb) > 1450* 10^(-9)
            sigma.epb(i)  = interp1(param_em(:,1), CrossLum, Lambda.pb(i) * 10^(9));            % ������� ������������� (��������� �������)
        else
            sigma.epb(i)  = 0;
        end
    end
end
    
