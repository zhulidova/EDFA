function [G, nf, P_out] = edfa_main(P_in, wl, n, low_gain_regime, L, T_c, r_edf, NA, splices) 
%% �������� ��������� �������
% ����������
% ph_const      - ��������� � ����������� �����������
% T             - ����������� � ���������
% splices.wdm_p - ������ �� wdm ��� ������� (datasheet), ���                         
% splices.wdm_s - ������ �� wdm ��� ������� (datasheet), ���                               
% splices.fiber - ������ �� ������ 
% n_sum         - ������������ � �^(-3)
% N             - ��������� � ��������� �������� ������� (N.S), ������� (N.PF � N.PB) � ASE (N.ASE)
% sigma         - ��������� � ��������� ��������� ������� (sigma.ES), ������� (sigma.EPF � sigma.EPB) � ASE (sigma.EASE)
%� ��������� ���������� ������� (sigma.AS), ������� (sigma.APF � sigma.APB) � ASE (sigma.AASE)
% psi           - ��������� � ��������� ������������ ������������� �� ������� ������� ��� ������� (psi.S),
%������� (psi.PF � psi.PB), ASE (psi.ASE) � �������������� ��������� ������������ (psi.NS, psi.PF, psi.PB, psi.ASE) 
% w_edf         - ��������� ������� �������� ��� ������� (w_edf.S), ������� (w_edf.PF � w_edf.PB) � ASE (w_edf.ASE) 

% �������
% concentration(n)              - ������� ��������� ������������ �� ppm � �^(-3)
% sigma_lum(T, wl, N, ph_const) - ������� ������� ������� ���������� � ������������� ��� ������ ����� ����� 
%��� ������������ �����������
% bessel(r_edf, wl, N, NA)      - ������� ������� ������� ������������ ������������� �� ������� ������� ���
%������ ���� ����� (������� �� ������� ���������� � �������� �������� ��������� �������)
% splice_loss(P_in.S,..., NA)   - ������� ������� �������� ������� ��������� ������� � ������� � ������ ������
%�� ������� � WDM
% odu2_unsaturated_gain_regime  - ������� ������� �������� ����������� ��� ���������� ��������� � �����������
%������� ������� 
% odu2                          - ������� ������� �������� ����������� ��� ���������� ��������� � ����� ������ 
% chord_method(L,..., ph_const) - ������� � �������� ������������� ��������� ������� ��� �������� ����������������
%������ ��� ������ �� ��������� �������� � ������ ���� 
% G(P_out, P_in.S)              - ������� ������� ������������ ��������
% nf(P_out, ph_const)           - ������� ������� ���-�������

ph_const.c             = 299792458;                                             % �������� �����    
ph_const.h             = 6.626E-34;                                             % ���������� ������
ph_const.tau           = 0.0102;                                                % ����� ����� �� ������� ������
ph_const.k             = 1.38 * 10^(-23);                                       % ���������� ���������

% ������ ���������� ��� ���
T                      = 273.15 + T_c;                                          % ����������� � ���������
n_sum                  = concentration(n);                                      % ������������ � �^(-3)
N.S                    = length(wl.S);                                          % ������ ������� ���� ���� �������
N.PF                   = length(wl.PF);                                         % ������ ������� ���� ���� �������� �������
N.ASE                  = length(wl.ASE);                                        % ������ ������� ���� ���� ASE
N.PB                   = length(wl.PB);                                         % ������ ������� ���� ��������� �������
sigma                  = sigma_lum(T, wl, ph_const);                            % ������ �������
[psi, w_edf]           = bessel(r_edf, wl, N, NA);                              % �������� ������� ������������� ������� ����

    %% ���� ���. ������  

[P_in.S, P_in.PF] = splice_loss(P_in.S, P_in.PF, wl, r_edf, splices, NA);      % �������� ������� �������� ������� � �������
                                                                              %� ������ ������ �� ������� � WDM

%% ������� �������� ���
    tStart = cputime;                                                           % ����� �������
    if isempty(P_in.PB) == 1 
% ������ ������ ������� ������� (������ ����)
        if low_gain_regime == 1
% ������� ���������� ������� ��������� (����� ������� �������)
            [z, P_out]     = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, wl, sigma, N, w_edf, n_sum, ph_const, P_in), [0: 0.1: L],[undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF)]);
        else
% ������� ��������� ������� ��������� (����� ������)
            [z, P_out]     = ode45(@(z,P) odu2(z, P, wl, sigma, psi, N, n_sum, ph_const, P_in), [0: 0.1: L],[undbm(P_in.S) undbm(P_in.ASEF) undbm(P_in.PF)]);
        end
        G                  = Gain(P_out(size(P_out,1), 1 : N.S)', P_in.S);      % ������ ��
        nf                 = NF(P_out(size(P_out,1), N.S + 1 : N.S + N.ASE),... % ������ ��
            G, wl, ph_const);                                                    
    else
% ������ ��������� ������� ��� ��������������� ������
        P_out              = chord_method(L, wl, sigma, psi, N, w_edf, n_sum, P_in, low_gain_regime, ph_const);
        G                  = Gain(P_out(size(P_out,1), 1 : N.S)', P_in.S);      % ������ ��
        nf                 = NF(P_out(size(P_out,1), N.S + 1 : N.S + N.ASE),... % ������ ��
            G, wl, ph_const);                                                   
    end
    P_out = P_out / undb(splices.fiber) / undb(splices.wdm_s);      % �������� ������� �������� ������� � �������
    time                   = cputime - tStart;                                  % ����� �������
    formatSpec             = '������� ��� ������ �� %f ������\n';                           
    fprintf(formatSpec,time)
end