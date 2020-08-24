%% ������� ������� ������� �������� ASE
%(Desurvire, appendix R)
% P_sat_P   - ������ �������� �������� ��������� ��� �������
% P_sat_ASE - ������ �������� �������� ��������� ��� ������� ASE
% n_sp      - ������ ������������� ���������� �������������, ����� ������ ������� ASE
% p_ase     - �������� ASE � ��
% P_ase     - �������� ASE � ���

function P_ase = ase_in(Lambda, ph_const, sigma, PinP, w_edf) 

P_sat_P   = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.pf.^2)...                  % �������� ��������� ��� �������
    ./ (Lambda.pf .* (sigma.apf + sigma.epf) * ph_const.tau);
P_sat_ASE = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.ase.^2)...                 % �������� ��������� ��� ASE
    ./ (Lambda.ase .* (sigma.aase + sigma.ease) * ph_const.tau);

n_sp      = 1 ./ (1 - sum(sigma.epf ./ sigma.apf) .* sigma.aase  ./ sigma.ease...     % ����������� ���������� �������������
    - sum((1 + sigma.epf./ sigma.apf) .* P_sat_P ./ PinP) .* sigma.aase  ./ sigma.ease);

if n_sp(1,1) < 0
    n_sp      = 1 ./ (1 - sum(sigma.epf ./ sigma.apf) .* sigma.aase  ./ sigma.ease...     % ����������� ���������� �������������
    - sum((1 + sigma.epf./ sigma.apf) .* P_sat_P ./ (3.5.*PinP)) .* sigma.aase  ./ sigma.ease);
end
    
P_ase     = 4 * n_sp .* ph_const.c * 0.1 * 10^(-9) ./ Lambda.ase.^2 .* P_sat_ASE * 10^(-16); % �������� ASE � ��

end