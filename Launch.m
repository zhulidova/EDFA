clear all;
% P_s_sum         - ��������� �������� �������, ���
% P_in.S          - ������ ������� ��������� ������� (��������� �������� / ���-�� ���� ����), ���
% P_in.PF         - ������ ��������� �������� ������� � 0, ���
% P_in.PB         - ������ ��������� ���������� ASE (�� �������� ������� ���� ����) � L, ���
%������ ������, ���� ������� ������ �������
% r_edf           - ������ ���������� ��������� �������, �
% NA              - �������� �������� ��������� �������    
% n               - ������������ ����� �����, ppm
% L               - ����� ��������� �������, �
% T_c             - ����������� �����, ��
% low_gain_regime - ������������ ����� ������ ������� ������� (���������� ������� ���������),
%0 - ����� ������ (��� �����������), 1 - ����� ������� �������

% �������
% edfa_main       - �������� ��������� ������� EDFA
% smooth_res      - ����������� �����������
% create_graf     - ���������� �������� 

P_s_sum = [-28];  % ������ ��������� ������� ��������� �������  
% ��������� ������� 
Lambda.s              = [1529.55 1532.70 1538.20 1543.80 1546.90 1551.70 1557.40 1560.60] * 10^(-9); % ����� ���� �������, ��
Pin.s                 = watt(P_s_sum) / length(Lambda.s) * ones(1,length(Lambda.s));        % ������� �������� �������, ���

% ��������� ASE 
Lambda_range          = [1528 1565] * 10^(-9);                           % ������ ASE, ��

% ��������� �������
Lambda.pf             = [1473 1480] * 10^(-9);                           % ����� ���� �������� �������, ��
Pin.pf                = watt([3.85 5.89]);

Lambda.pb             = [] * 10^(-9);                                    % ����� ���� ��������� �������, ��
                                                                         % (������ ������, ���� ������� ������ ��������)
Pin.pb                = [];                                              % ������� �������� ��������� �������, ���

L    = 9;
N0.r                  = 25;
N0.ase                = 20;
N0.z                  = 50;
type                  = 'OFS980';

T_c = 25;

%% �������� ��������� �������
[z, P, Gain, OSNR, NF, SPase]   = edfa(Pin, Lambda, Lambda_range, L, N0, type);
% ����������� ����������� �������������
[wl_smooth.Model, G_smooth.Model, nf_smooth.Model] = smooth_res(Lambda.s', Gain, NF);   

% ������ ����������������� ������
EXP   = exp_res(P_s_sum, T_c, dbm(sum(Pin.pf)));

[wl_smooth.Exp, G_smooth.Exp, nf_smooth.Exp]  = smooth_res(EXP(1,:), EXP(2,:), EXP(3,:));

%% ���������� ��������
create_graf(wl_smooth.Model*10^9, G_smooth.Model, "model", 'Gain', P_s_sum, Pin, 1, 1); 
hold on;
create_graf(wl_smooth.Exp, G_smooth.Exp, "exp", 'Gain', P_s_sum, Pin, 1, 1); 
hold on;
create_graf(wl_smooth.Model*10^9, nf_smooth.Model, "model", 'NF', P_s_sum, Pin, 2, 1); 
hold on;
create_graf(wl_smooth.Exp, nf_smooth.Exp, "exp", 'NF', P_s_sum, Pin, 2, 1);
hold on;
