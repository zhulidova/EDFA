%% ������� ���������
% �	Pin � ���������, ��������� �� �������� ������� ��������� ��������, �������� � ��������� �������.
% �	Lambda � ���������, ��������� �� �������� �������� ��������������� ���� ���� ��� ��������, �������� � ��������� �������. 
% �	Lambda_range � ������, �������� �������� ������ ������� ��������� ���� ���� ��� ������� ����   
% �	L � ����� ��������� ������� 
% �	Type � ��������, ����������� �� ������������ ��� ������� (��������� ��������� �������). 
% �	N0 � ���������, ���������� �������� ���������� ����� ��������� ��������� ���� ���� ��� ������� ����, ���������� ����� ������������� 
% ������� ���������� ������� (��� �������������� �� �������), ���������� ����� ������������� ����� ������� (��������������� ��� �������������� �� 
% ����� �������). 
%
% ������ 
% Pin.s = []*10^(-3);     % �������� ��������, ��
% Lambda.s = []*10^(-9);  % ����� ���� ��������, �
% Pin.pf = []*10^(-3);    % �������� �������� �������, ��
% Lambda.pf = []*10^(-9); % ����� ���� �������� �������, �
% Pin.pb = []*10^(-3);    % �������� ��������� �������, ��
% Lambda.pb = []*10^(-9); % ����� ���� ��������� �������, �
% Lambda_range = [1528 1565]*10^(-9); % �������� ���� ���� ��� ������� ����, �
% L = [];                 % ����� ��������� �������, �
% N.ase = 20;             % ����� ����� �� ����� ����� (��� ������� ����) 
% N.z = 50;               % ����� ����� �� z (��� ��������������)
% N.r = 25;               % ����� ����� �� r
% type = 'OFS980';        % ��� ��������� �������

%% �������� ���������
% �	z � ������ �������� ����� �������, �� ������� ������������ ������ ����������. 
% �	P � ���������, ��������� �� �������� �������� ��������� ��������, �������� � ��������� �������,
% �������� � ��������� �����, ��������� ������������� ��������� �� z. �������� � ��. 
% �	Gain � ������ �������� ������������� �������� ��� ���������������
% ������� ��������, dB
% �	OSNR � ������ �������� ����������� ������/��� ��� ��������������� ������� ��������, dB  
% �	NF � ������ �������� ���-������� ��� ��������������� ������� ��������.
% dB
% �	SPase � ������ ���� (����������� �������� ���� �� ����� �����)

%% �������
% concentration(...) - ������� ��������� ������������ �� ppm � �^(-3)
% sigma_lum(...)     - ������� ������� ������� ���������� � ������������� ��� ������ ����� ����� 
%��� ������������ �����������
% bessel(...)        - ������� ������� ������� ������������ ������������� �� ������� ������� ���
%������ ���� ����� (������� �� ������� ���������� � �������� �������� ��������� �������)
% splice_loss(...)   - ������� ������� �������� ������� ��������� ������� � ������� � ������ ������
%�� ������� � WDM
% ase_in(...)        - ������� ������� ��������� �������� �������� ASE
% odu2_unsaturated_gain_regime(...) - ������� ������� �������� ����������� ��� ���������� ��������� � �����������
%������� ������� 
% odu2(...)          - ������� ������� �������� ����������� ��� ���������� ��������� � ����� ������ 
% chord_method(...)  - ������� � �������� ������������� ��������� ������� ��� �������� ����������������
%������ ��� ������ �� ��������� �������� � ������ ���� 
% gain(...)          - ������� ������� ������������ ��������
% nf(...)            - ������� ������� ���-�������

function [z, P, Gain, OSNR, NF, SPase] = edfa(Pin, Lambda, Lambda_range, L, N0, type) 

% ���������� ���������
ph_const.c      = 299792458;                                         % �������� �����    
ph_const.h      = 6.626E-34;                                         % ���������� ������
ph_const.tau    = 0.0102;                                            % ����� ����� �� ������� ������
ph_const.k      = 1.38 * 10^(-23);                                   % ���������� ���������
ph_const.eV     = 1.602 * 10^(-19);                                  % 1 �� � ��

% ��������� ��������� �������
if strcmp(type,'OFS980')==1
    r_edf       = 2.9E-6 / 2;                                        % ������ ����������, �
    NA          = 0.26;                                              % �������� �������� ��������� �������
    n           = 70;                                                % ������������ ����� �����, ppm
end 

% ��������� �������
T_c             = 25;                                                 % ����������� �����, ��
splices.wdm_p   = 2.35;                                              % ������ �� wdm ��� ������� (datasheet), ���                         
splices.wdm_s   = 0.2;                                               % ������ �� wdm ��� ������� (datasheet), ���                  
splices.fiber   = 0.5;                                               % ������ �� ������
low_gain_regime = 0;                                                 % ������������ ����� ������ ������� �������
% 1 - ����� ������� ������� % 0 - ����� ������

% ������ ���������� ��� ���
T               = 273.15 + T_c;                                      % ����������� � ���������
n_sum           = concentration(n);                                  % ������������ � �^(-3)
N.s             = length(Lambda.s);                                  % ������ ������� ���� ���� �������
N.pf            = length(Lambda.pf);                                 % ������ ������� ���� ���� �������� �������
N.ase           = N0.ase;                                            % ������ ������� ���� ���� ASE
N.pb            = length(Lambda.pb);                                 % ������ ������� ���� ��������� �������
Lambda.ase      = ase_diskr(N.ase,Lambda_range);                     % ������ ��������� �������� �������� ASE
sigma           = sigma_lum(T, Lambda, ph_const);                    % ������ �������
[psi, w_edf]    = bessel(r_edf, Lambda, N, NA, N0.r);                % �������� ������� ������������� ������� ����

%% ���� ���. ������  

[Pin.s, Pin.pf] = splice_loss(dbm(Pin.s), dbm(Pin.pf), Lambda, r_edf, splices, NA); 
Pin.asef        = ase_in(Lambda, ph_const, sigma, Pin.pf, w_edf);

%% ������� �������� ���
    tStart = cputime;                   % ����� �������
    if isempty(Pin.pb) == 1 
    % ������ ������ ������� ������� (������ ����)
        if low_gain_regime == 1
        % ������� ���������� ������� ��������� (����� ������� �������)
            [z, Pout]   = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, Lambda, sigma, N, w_edf, n_sum, ph_const, Pin), [0: L/N0.z: L],[Pin.s Pin.asef  Pin.pf]);
        else
        % ������� ��������� ������� ��������� (����� ������)
            [z, Pout]   = ode45(@(z,P) odu2(z, P, Lambda, sigma, psi, N, n_sum, ph_const, Pin, r_edf, N0), [0: L/N0.z: L],[Pin.s Pin.asef Pin.pf]);
        end
    else
        % ������ ��������� ������� ��� ��������������� ������
        Pin.aseb        = ase_in(Lambda,ph_const, sigma,Pin.pb, w_edf);
        [z, Pout]       = chord_method(L, Lambda, sigma, psi, N, w_edf, n_sum, Pin, low_gain_regime, ph_const, r_edf, N0); 
    end
    Pout       = Pout / undb(splices.fiber); % �������� �������� �������� ������� � �������
    
    z          = z';   % ������ �������� z
    % P - ��������� ������������� ��������� �� z
    P.s        = Pout(:, 1: N.s)';                      % ������������� ��������� ��������
    Pspec.s    = P.s(:, size(P.s,2));                   % ������ �������� ������� � ����� �������
    P.asef     = Pout(:, N.s+1: N.s+N.ase)';            % ������������� �������� �������� �����
    Pspec.asef = P.asef(:, size(P.asef,2));             % ������ �������� �������� �����
    P.pf       = Pout(:, N.s+N.ase+1: N.s+N.ase+N.pf)'; % ������������� ��������� �������� �������
    Pspec.pf   = P.pf(:, size(P.pf,2));                 % ������ �������� �������� �������

    if isempty(Pin.pb) == 0
        P.aseb     = Pout(:, N.s+N.ase+N.pf+1: N.s+2*N.ase+N.pf)';        % ������������� �������� ��������� �����
        Pspec.aseb = P.aseb(:, size(P.aseb,2));                           % ������ �������� ��������� �����
        P.pb       = Pout(:, N.s+2*N.ase+N.pf+1: N.s+2*N.ase+N.pf+N.pb)'; % ������������� ��������� ��������� �������
        Pspec.pb   = P.pb(:, size(P.pb,2));                               % ������ �������� ��������� �������

    end
    
    Gain         = gain(Pspec.s, Pin.s');                           % ������ Gain
    [OSNR, NF]   = nf(Pspec.asef, Gain, Lambda, ph_const, Pspec.s); % ������ NF � OSNR_out
    SPase.lambda = Lambda.ase';                                     % ������ ASE
    SPase.p      = dbm(Pspec.asef);
    
    time         = cputime - tStart;                                % ����� �������
    formatSpec   = '������� ��� ��� EDFA ������ �� %f ������\n';                           
    fprintf(formatSpec,time)
    
end