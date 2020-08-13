clear all; 
%% �������� ���� � �������� �������
% wl.S            - ������ ���� ���� �������, ��
% wl.ASE          - ������ ���� ���� ASE, ��
% wl.PF           - ������ ���� ���� �������� �������, ��
% P_s_sum         - ��������� �������� �������, ���
% P_in.S          - ������ ������� ��������� ������� (��������� �������� / ���-�� ���� ����), ���
% P_in.PF         - ������ ��������� �������� ������� � 0, ���
% P_in.ASEF       - ������ ��������� ��������� ASE (�� �������� ������� ���� ����) � 0, ��� 
% P_in.ASEB       - ������ ��������� ��������� ������� � L, ���
% P_in.PB         - ������ ��������� ���������� ASE (�� �������� ������� ���� ����) � L, ���
%������ ������, ���� ������� ������ �������
% r_edf           - ������ ���������� ��������� �������, �
% NA              - �������� �������� ��������� �������    
% n               - ������������ ����� �����, ppm
% L               - ����� ��������� �������, �
% T_c             - ����������� �����, ��
% low_gain_regime - ������������ ����� ������ ������� ������� (���������� ������� ���������),
%0 - ����� ������ (��� �����������), 1 - ����� ������� �������

P_s_sum = 2;    % ������ ��������� ������� ��������� �������
                                
% ��������� �������
wl.S                   = [1529.55 1532.52 1537.24 1541.96 1547.56 1552.36 1556.36 1560.44];% ����� ���� �������, ��
P_in.S                 = dbm(undbm(P_s_sum) / length(wl.S)) * ones(1,length(wl.S));        % ������� �������� �������, ���

% ��������� ASE 
wl.ASE                 = [1528 : 0.1 : 1565];                                              % ������ ASE, ��

% ��������� �������
wl.PF                  = [1473 1480];                                                      % ����� ���� �������� �������, ��
P_in.PF                = [13.85 15.89];                                                    % ������� �������� �������� �������, ���
wl.PB                  = [1473 1480];                                                      % ����� ���� ��������� �������, ��
% (������ ������, ���� ������� ������ ��������)
P_in.PB                = [];                                                               % ������� �������� ��������� �������, ���

% ��������� ��������� ������� � �������
r_edf                  = 2.9E-6 / 2;                                                       % ������ ����������, �
NA                     = 0.26;                                                             % �������� �������� ��������� �������
n                      = 140;                                                              % ������������ ����� �����, ppm
L                      = 9;                                                                % ����� ��������� �������, �
T_c                    = 100;                                                               % ����������� �����, ��
splices.wdm_p          = 0.15;                                                             % ������ �� wdm ��� ������� (datasheet), ���                         
splices.wdm_s          = 0.2;                                                              % ������ �� wdm ��� ������� (datasheet), ���                  
splices.fiber          = 0.5;                                                              % ������ �� ������
low_gain_regime        = 0;                                                                % ������������ ����� ������ ������� �������
                                                                                           % 1 - ����� ������� ������� % 0 - ����� ������
%% �������� ��������� �������
[G, nf, Pout]          = edfa_main(P_in, wl, n, low_gain_regime, L, T_c, r_edf, NA, splices);

% ����������� ����������� �������������
[wl_smooth.Model, G_smooth.Model, nf_smooth.Model] = smooth_res(wl.S, G, nf);   

% ������ ����������������� ������
EXP                    = exp_res(P_s_sum, T_c, dbm(sum(undbm(P_in.PF))));

[wl_smooth.Exp, G_smooth.Exp, nf_smooth.Exp]       = smooth_res(EXP(1,:), EXP(2,:), EXP(3,:));

delta_G = G - EXP(2,:);
delta_nf = nf - EXP(3,:);
%% ���������� ��������
create_graf(wl_smooth.Model, G_smooth.Model, "model", 'Gain', P_s_sum, 1, 1); 
hold on;
create_graf(wl_smooth.Exp, G_smooth.Exp, "exp", 'Gain', P_s_sum, 1, 1); 
hold on;
create_graf(wl_smooth.Model, nf_smooth.Model, "model", 'NF', P_s_sum, 2, 1); 
hold on;
create_graf(wl_smooth.Exp, nf_smooth.Exp, "exp", 'NF', P_s_sum, 2, 1); 
