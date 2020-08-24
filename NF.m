%% ������ ������� ���-�������
% ������ �� Giles, Desurvire, 1991, Propagation of signal and noise in
% concatenated erbium-doped fiber optical amplifiers

% ������� ���-������� (Desurvire):
% (1 + 2 .* n_sp .* (undb(G) - 1)) ./ undb(G)            - �������� ���
% M .* (n_sp .* (undb(G) - 1)).^2 ./(undb(G)).^2/n_mean  - ���������-���������� ��� ������
% M .* n_sp .* (undb(G) - 1)./(undb(G)).^2/n_mean        - ������-���������� ��� ������

function [OSNR, NF] = nf(P_out_ase,G, Lambda, ph_const, P_out_s)

P_ase  = interp1(Lambda.ase,P_out_ase',Lambda.s);          % ���������� �������� ASE �� ����� ����� �������

% M      = 1;                                              % ����� ���
% n_mean = 2;                                              % ������� ����� ������� ��� z = 0, ������������� ��������

n_sp   = P_ase' ./ (2 * (undb(G) - 1) * ph_const.h...      % ����������� ���������� �������������
    * ph_const.c^2 * 0.1 * 10^(-9) ./ Lambda.s'.^3);   
f      = (1 + 2 .* n_sp .* (undb(G) - 1)) ./ undb(G);      % ���-������ � �������� ��������
NF     = db(f);                                            % ���-������ � ��
OSNR   = db(P_out_s ./ P_ase');
end