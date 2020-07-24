%% ������� ������� �������� ������� �������� ������� � ������� � ������ ������ �� ������� � WDM
function [P_s, P_p] = splice_loss(P_s, P_p, wl, r_edf, splices, NA_edf)
    

% ��������� ������� SMF-28E
r_smf   = 8.20E-6/2;                                 % ������ ���������
NA_smf  = 0.14;                                      % �������� ��������
V_p_smf = 2 * pi * r_smf * NA_smf * 10^9 ./ wl.PF;   % ������������� ������� (�������)
V_s_smf = 2 * pi * r_smf * NA_smf * 10^9 ./ wl.S;    % ������������� ������� (������)
 
% ������� ������� SMF-28E
w_p_smf = r_smf .* (0.65 + (1.619 ./ V_p_smf.^1.5) + (2.879 ./ V_p_smf.^6)); 
w_s_smf = r_smf .* (0.65 + (1.619 ./ V_s_smf.^1.5) + (2.879 ./ V_s_smf.^6)); 
 
% ��������� ��������� �������
V_p_edf  =  2 * pi * r_edf * NA_edf * 10^9 ./ wl.PF; % ������������� ������� (�������)
V_s_edf  =  2 * pi * r_edf * NA_edf * 10^9 ./ wl.S;  % ������������� ������� (������)
 
% ������� ������� ��������� �������
w_p_edf  =  r_edf .* (0.65 + (1.619 ./ V_p_edf.^1.5) + (2.879 ./ V_p_edf.^6));
w_s_edf  =  r_edf .* (0.65 + (1.619 ./ V_s_edf.^1.5) + (2.879 ./ V_s_edf.^6)); 

% L_mismatch_p =   db(4 ./ ((w_p_smf ./ w_p_edf).^2 + (w_p_edf ./...
%  w_p_smf).^2));
% L_mismatch_s =   db(4 ./ ((w_s_smf ./ w_s_edf).^2 + (w_s_edf ./...
%  w_s_smf).^2)); 

% �������� �������� �������� �� ����� � �������� �������
P_s  = P_s - splices.fiber;           % �������� �������� �������
P_p  = P_p - splices.wdm_p - splices.fiber;           % �������� �������� �������
end