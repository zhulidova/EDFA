%% ������� ������� ������ � ������� ������ ���� (�������, ����� �������)
% � ������� ������ ���� ��������� �������� ��������� ��������� ������� � ���������� ���� � ������ �����
%(�������� ������ � ���������� ��������� � ������ ����)���������� ������ - ������ �� ������ EDFA/��������. ����� �� 27.06)
function [z, P_out] = chord_method(L, Lambda, sigma, psi, N,w_edf, n_sum, Pin, low_gain_regime, ph_const, r_edf, N0)
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
epsilon = 1e-9;
        ksi1(1,:)         = zeros(1, N.ase);            % ��������� ������� ��� ������ ���� (ASE,������ ����� ��������)
        ksi2(1,:)         = ones(1, N.ase) * 10^(-8);   % ��������� ������� ��� ������ ���� (ASE,����� ����� ��������)
        lambda1(1,:)      = zeros(1, N.pb)* 0.001;      % ��������� ������� ��� ������ ���� (P_p_backward,������ ����� ��������)
        lambda2(1,:)      = ones(1, N.pb) * 0.005;      % ��������� ������� ��� ������ ���� (P_p_backward,����� ����� ��������)
        lambda            = ones(1, N.pb) * 1;
%% ���� ��� ������ ������������ ������� ��� ������ ����
     while abs(sum(lambda1) - sum(lambda)) > epsilon
         if min(lambda1) < 0
             lambda1(1,:) = ones(1, N.pb) * 3*10^(-6);  % ��������� ������� ��� ������ ���� (P_p_backward,������ ����� ��������)
             lambda2(1,:) = ones(1, N.pb) * 2*10^(-5);
         end
         if low_gain_regime == 1 
            % ������� ���������� ������� ��������� (����� ������� �������) 
            [z, P_out1]   = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, Lambda, sigma, N,w_edf, n_sum, ph_const, Pin), [0: L/N0.z: L], [Pin.s Pin.asef Pin.pf ksi1 lambda1]);  
            [z, P_out2]   = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, Lambda, sigma, N,w_edf, n_sum, ph_const, Pin), [0: L/N0.z: L], [Pin.s Pin.asef Pin.pf ksi2 lambda2]); 
         else 
            % ������� ��������� ������� ��������� (����� ������)
            [z, P_out1]   = ode45(@(z,P) odu2(z, P, Lambda, sigma, psi, N, n_sum, ph_const, r_edf, N0), [0: L/N0.z: L],[Pin.s Pin.asef Pin.pf ksi1 lambda1]);  
            [z, P_out2]   = ode45(@(z,P) odu2(z, P, Lambda, sigma, psi, N, n_sum, ph_const, r_edf, N0), [0: L/N0.z: L],[Pin.s Pin.asef Pin.pf ksi2 lambda2]); 
         end
        
        % ������ �������� ��������� ������� � ���������� ASE � ����� �����
        P_ASEB1   = P_out1(size(P_out1,1), N.s+N.ase+N.pf+1: N.s+N.ase+N.pf+N.ase);                       % ���������� �������� ASE � ����� �������, 1
        P_ASEB2   = P_out2(size(P_out2,1), N.s+N.ase+N.pf+1: N.s+N.ase+N.pf+N.ase);                       % ���������� �������� ASE � ����� �������, 2
        P_PB1     = P_out1(size(P_out1,1), N.s+2*N.ase+N.pf+1: N.s+N.ase+N.pf+N.ase+N.pb);                % ���������� �������� P_p_backward � ����� �������, 1
        P_PB2     = P_out2(size(P_out2,1), N.s+2*N.ase+N.pf+1: N.s+N.ase+N.pf+N.ase+N.pb);                % ���������� �������� P_p_backward � ����� �������, 2
        
        % ������������ ����������� ������� ��� ������ ����
        ksi       = ksi1 + (ksi2 - ksi1) .* (Pin.aseb - P_ASEB1) ./ (P_ASEB2 - P_ASEB1);
        lambda    = lambda1 + (lambda2 - lambda1) .* (Pin.pb - P_PB1) ./ (P_PB2 - P_PB1);
        m         = lambda1;
        lambda1   = lambda;
        lambda    = m;
        ksi1      = ksi;    
     end
     
     % ������� ���������� ������ ����
     if low_gain_regime == 1
        % ������� ���������� ������� ��������� (����� ������� �������)
        [z, P_out]  = ode45(@(z,P) odu2_unsaturated_gain_regime(z, P, Lambda, sigma, N,w_edf, n_sum, ph_const, Pin), [0: L/N0.z: L],[Pin.s Pin.asef Pin.pf ksi1 lambda1]);
        
     else
         % ������� ��������� ������� ��������� (����� ������)
        [z, P_out]  = ode45(@(z,P) odu2(z, P, Lambda, sigma, psi, N, n_sum, ph_const, r_edf, N0), [0: L/N0.z: L], [Pin.s Pin.asef Pin.pf ksi1 lambda1]);
     end
end