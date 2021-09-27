clear all; clc;
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
P_p(1,:) = [1.85 3.89];
P_p(2,:) = [2.85 4.89];
P_p(3,:) = [3.85 5.89];
P_p(4,:) = [5.85 7.89];
P_p(5,:) = [7.85 9.89];
P_p(6,:) = [9.85 11.89];
P_p(7,:) = [11.85 13.89];
P_p(8,:) = [13.85 15.89];

Len = [6,7,8,9,10,11,12,13,14];
Temp = [-60, -50, -40, -20, 0, 25, 40, 70];
P_s = [-8 -3 2 7 10];
for k = 1:size(Temp,2)
for j = 1:size(P_p,1)
for i = 1:size(P_s,2)
    
P_s_sum = P_s(i);  % ������ ��������� ������� ��������� �������  
% ��������� ������� 
Lambda.s              = [1550.1	1551.7	1554.1	1557.33	1559.0	1560.58	1562.2	1564.65]* 10^(-9); % ����� ���� �������, ��
Pin.s                 = watt(P_s_sum) / length(Lambda.s) * ones(1,length(Lambda.s));        % ������� �������� �������, ���

% ��������� ASE 
Lambda_range          = [1528 1565] * 10^(-9);                           % ������ ASE, ��

% ��������� �������
Lambda.pf             = [1473 1480] * 10^(-9);                           % ����� ���� �������� �������, ��
Pin.pf                = watt(P_p(j,:));

Lambda.pb             = [];                                               % ����� ���� ��������� �������, ��
                                                                         % (������ ������, ���� ������� ������ ��������)
Pin.pb                = [];                                              % ������� �������� ��������� �������, ���

L    = 9;
N0.r                  = 25;
N0.ase                = 50;
N0.z                  = 50;
type                  = 'Pnppk_a2';

T_c = Temp(k);

%% �������� ��������� �������
% [z, P, Gain, OSNR, NF, SPase]   = edfa(Pin, Lambda, Lambda_range, L, N0, type, T_c);
% Pout = dbm(P.s);


% ����������� ����������� �������������
% [wl_smooth.Model, G_smooth.Model, nf_smooth.Model] = smooth_res(Lambda.s', Gain, NF); 

% ������ ����������������� ������
EXP   = exp_res(P_s_sum, T_c, dbm(sum(Pin.pf)),L);
  
[wl_smooth.Exp, G_smooth.Exp, nf_smooth.Exp]  = smooth_res(EXP(1,:), EXP(2,:), EXP(3,:));


%% ���������� ��������
% create_graf(Lambda.s*10^9, Gain, "model", 'Gain', P_s_sum, Pin, 1, i); 
% hold on;
% create_graf(EXP(1,:), EXP(2,:), "exp", 'Gain', P_s_sum, Pin, 1, i, T_c); 
% hold on;
% create_graf(Lambda.s*10^9, NF, "model", 'NF', P_s_sum, Pin, 2, i); 
% hold on;
create_graf(wl_smooth.Exp, nf_smooth.Exp, "exp", 'NF', P_s_sum, Pin, 2, i, T_c);
hold on;
end
% Filename1 = "F_G" + sprintf("%.0f", T_c)+"_"+sprintf("%.0f", dbm(sum(Pin.pf)))+".emf";
Filename2 = "F_NF" + sprintf("%.0f", T_c)+"_"+sprintf("%.0f", dbm(sum(Pin.pf)))+".emf";
% saveas(gcf,Filename1);
saveas(gcf,Filename2);
hold off;
end
end
% plot(Len,G);
% hold on;
% plot(Len,G_exp);
% semilogx(watt(P_p)*10^3,T);
% hold on;
% P_in_exp = [-9.96 -7.99 -5.03 -2.97 0.04 2.02 4.97 6.99];
% P_out_exp = [-11.4400000000000 -9.28000000000000 -6.10000000000000 -3.85000000000000 -0.630000000000000 1.44000000000000 4.46000000000000 6.56000000000000];
% T_exp = watt(P_out_exp+0.33)./watt(P_in_exp);
% semilogx(watt(P_in_exp)*10^3,T_exp);
