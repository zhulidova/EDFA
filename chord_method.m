%% ������� ������� ������ � ������� ������ ���� (�������, ����� �������)
% � ������� ������ ���� ��������� �������� ��������� ��������� ������� � ���������� ���� � ������ �����
%(�������� ������ � ���������� ��������� � ������ ����)���������� ������ - ������ �� ������ EDFA/��������. ����� �� 27.06)
function P_out = chord_method(L, wl, sigma, psi, N,w_edf, n_sum, P_in, low_gain_regime, ph_const, r_edf)
% ���������� 
% epsilon           - �������� �������� ������������� ������
% ksi1, ksi2        - �������� ��������� ������� ��� ASE
% lamda1, lamda2    - �������� ��������� ������� ��� �������
% lamda             - ���������� �������� lamda1 
% P_ASEB1 � P_ASEB2 - �������� �������� ��������� ����� � ����� �����
% P_PB1 � P_PB2     - �������� �������� ��������� ������� � ����� �����

% �������
% odu2_unsaturated_gain_regime  - ������� ������� �������� ����������� ��� ���������� ��������� � �����������
%������� ������� 
% odu2                          - ������� ������� �������� ����������� ��� ���������� ��������� � ����� ������
epsilon = 1e-6;
        ksi1(1,:)           = zeros(1, N.ASE);            % ��������� ������� ��� ������ ���� (ASE,������ ����� ��������)
        ksi2(1,:)           = ones(1, N.ASE) * 10^(-9);   % ��������� ������� ��� ������ ���� (ASE,����� ����� ��������)
        lamda1(1,:)         = zeros(1, N.PB);             % ��������� ������� ��� ������ ���� (P_p_backward,������ ����� ��������)
        lamda2(1,:)         = ones(1, N.PB) * 0.0005;       % ��������� ������� ��� ������ ���� (P_p_backward,����� ����� ��������)
        lamda               = ones(1, N.PB) * 1;
%% ���� ��� ������ ������������ ������� ��� ������ ����
     while abs(sum(lamda1) - sum(lamda)) > epsilon
         
         if low_gain_regime == 1 
            % ������� ���������� ������� ��������� (����� ������� �������) 
            [z, P_out1]     = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, wl, sigma, N,w_edf, n_sum, ph_const, P_in), [0: 0.1: L], [undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF) ksi1 lamda1]);  
            [z, P_out2]     = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, wl, sigma, N,w_edf, n_sum, ph_const, P_in), [0: 0.1: L], [undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF) ksi2 lamda2]); 
         else 
            % ������� ��������� ������� ��������� (����� ������)
            [z, P_out1]     = ode45(@(z,P) odu2(z, P, wl, sigma, psi, N, n_sum, ph_const, P_in, r_edf), [0: 0.1: L],[undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF) ksi1 lamda1]);  
            [z, P_out2]     = ode45(@(z,P) odu2(z, P, wl, sigma, psi, N, n_sum, ph_const, P_in, r_edf), [0: 0.1: L],[undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF) ksi2 lamda2]); 
         end
        
        % ������ �������� ��������� ������� � ���������� ASE � ����� �����
        P_ASEB1     = P_out1(size(P_out1,1), N.S + N.ASE + N.PF + 1 : N.S + N.ASE + N.PF + N.ASE);                           % ���������� �������� ASE � ����� �������, 1
        P_ASEB2     = P_out2(size(P_out2,1), N.S + N.ASE + N.PF + 1 : N.S + N.ASE + N.PF + N.ASE);                           % ���������� �������� ASE � ����� �������, 2
        P_PB1       = P_out1(size(P_out1,1), N.S + 2 * N.ASE + N.PF + 1 : N.S + N.ASE + N.PF + N.ASE + N.PB);                % ���������� �������� P_p_backward � ����� �������, 1
        P_PB2       = P_out2(size(P_out2,1), N.S + 2 * N.ASE + N.PF + 1 : N.S + N.ASE + N.PF + N.ASE + N.PB);                % ���������� �������� P_p_backward � ����� �������, 2
        
        % ������������ ����������� ������� ��� ������ ����
        ksi         = ksi1 + (ksi2 - ksi1) .* (undbm(P_in.ASEB) - P_ASEB1) ./ (P_ASEB2 - P_ASEB1);
        lamda       = lamda1 + (lamda2 - lamda1) .* (undbm(P_in.PB) - P_PB1) ./ (P_PB2 - P_PB1);
        m           = lamda1;
        lamda1      = lamda;
        lamda       = m;
        ksi1        = ksi;    
     end
     
     % ������� ���������� ������ ����
     if low_gain_regime == 1
        % ������� ���������� ������� ��������� (����� ������� �������)
        [z, P_out]  = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, wl, sigma, N,w_edf, n_sum, ph_const, P_in), [0: 0.1: L],[undbm(P_in.S)...
            undbm(P_in.ASEF) undbm(P_in.PF) ksi1 lamda1]);
        
     else
         % ������� ��������� ������� ��������� (����� ������)
        [z, P_out]  = ode45(@(z,P) odu2(z, P, wl, sigma, psi, N, n_sum, ph_const, P_in, r_edf), [0: 0.1: L], [undbm(P_in.S) undbm(P_in.ASEF)...
            undbm(P_in.PF) ksi1 lamda1]);
     end
end