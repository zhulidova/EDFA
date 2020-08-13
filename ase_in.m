%% ������� ������� ������� �������� ASE
%(Desurvire, appendix R)
% P_sat_P   - ������ �������� �������� ��������� ��� �������
% P_sat_ASE - ������ �������� �������� ��������� ��� ������� ASE
% n_sp      - ������ ������������� ���������� �������������, ����� ������ ������� ASE
% p_ase     - �������� ASE � ��
% P_ase     - �������� ASE � ���

function P_ase = ase_in(wl,ph_const, sigma,P_in, w_edf) 

P_sat_P   = (ph_const.h * ph_const.c * 10^34 * pi .* w_edf.PF.^2)...               % �������� ��������� ��� �������
    ./ (wl.PF .* (sigma.APF + sigma.EPF) * ph_const.tau);
P_sat_ASE = (ph_const.h * ph_const.c * 10^34 * pi .* w_edf.ASE.^2)...              % �������� ��������� ��� ASE
    ./ (wl.ASE .* (sigma.AASE + sigma.EASE) * ph_const.tau);

n_sp      = 1 ./ (1 - sum(sigma.EPF./ sigma.APF) .* sigma.AASE  ./ sigma.EASE...   % ����������� ���������� �������������
    - sum((1 +sigma.EPF./ sigma.APF) .*P_sat_P./P_in).* sigma.AASE  ./ sigma.EASE);

p_ase     = 4 * n_sp .* ph_const.c * 0.1 * 10^9 ./ wl.ASE.^2 .* P_sat_ASE * 10^(-16); % �������� ASE � ��
P_ase     = dbm(p_ase);        % �������� ASE � ���
end