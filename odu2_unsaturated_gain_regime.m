function f = odu2_unsaturated_gain_regime(z, P, wl, sigma, N,w_edf, n_sum, ph_const, P_in)
%% ���������� ������� ��� (����� ������� �������)
% ������ ���������� ��� ������ �������� ������� ������������� ������������ (����� -
%������ �� ������ EDFA/��������. ����� �� 02.06)

% ����������
% P_sat   - ��������� ��������� ��������� ��� ������� (P_sat.S), ������� (P_sat.PF � P_sat.PB) � ASE (P_sat.ASE)
% p       - ��������� ����������� ��������� ������� (p.S), ������� (p.PF � p.PB) � ASE (p.ASEF � p.ASEB)
% overlap - ��������� ���������� ����������: �������� (overlap.ES, overlap.EASE), ���������� (overlap.AS, overlap.AASE)
%overlap.sum_X - ���������� ��� ������� ���������� ����������
% g       - ��������� ���������������� ������������� �������� (g.ES � g.EASE) � ���������� (g.AS � g.AASE)
% alpha   - ��������� ������������ ��������� �������� (alpha.PF) � ��������� (alpha.PB) �������
% f       - �������� ����������� � ���������� ����������

P_sat.S             = (ph_const.h * ph_const.c * 10^34 * pi .* w_edf.S.^2)...   % �������� ��������� (������)
    ./ (wl.S .* (sigma.AS + sigma.ES) * ph_const.tau);
P_sat.PF            = (ph_const.h * ph_const.c * 10^34 * pi .* w_edf.PF.^2)...  % �������� ��������� (�������� �������)
    ./ (wl.PF .* (sigma.APF + sigma.EPF) * ph_const.tau);                                   
P_sat.PB            = (ph_const.h * ph_const.c * 10^34 * pi .* w_edf.PB.^2)...  % �������� ��������� (��������� �������)
    ./ (wl.PB .* (sigma.APB + sigma.EPB) * ph_const.tau);                                   
P_sat.ASE           = (ph_const.h * ph_const.c * 10^34 * pi .* w_edf.ASE.^2)... % �������� ��������� (ASE)
    ./ (wl.ASE .* (sigma.AASE + sigma.EASE) * ph_const.tau);

p.S                 = P(1: N.S, 1) ./ P_sat.S';                                 % ����������� �������� �������
p.ASEF              = P(N.S+1:N.S+N.ASE, 1) ./ P_sat.ASE';                      % ����������� �������� ��������� ASE
p.PF                = P(N.S+N.ASE+1:N.S+N.ASE+N.PF,1) ./ P_sat.PF';             % ����������� �������� �������� �������
overlap.sum_ES      = sum(p.S' .* sigma.ES ./ (sigma.AS + sigma.ES)) +...       % ��������� � ��������� ��� ������� � ASE
    sum(p.ASEF' .* sigma.EASE ./ (sigma.AASE + sigma.EASE));
overlap.sum_all     = sum(p.S) + sum(p.ASEF) + sum(p.PF);                       % ����� ���� ����������� ���������
overlap.sum_EP      = sum(p.PF' .* sigma.EPF ./ (sigma.APF + sigma.EPF));      % ��������� � ��������� ��� �������
overlap.sum_AP      = sum(p.PF' .* sigma.APF ./ (sigma.APF + sigma.EPF));       % ��������� � ����������� ��� �������

if isempty(P_in.PB) == 0
    p.ASEB          = P(N.S+N.ASE+N.PF+1:N.S+N.ASE+N.PF+N.ASE) ./ P_sat.ASE';   % ����������� �������� ��������� ASE
    p.PB            = P(N.S+N.ASE+N.PF+N.ASE+1:N.S+N.ASE+N.PF+N.ASE+N.PB) ...   % ����������� �������� ��������� �������
        ./ P_sat.PB';
    overlap.sum_ES  = overlap.sum_ES + sum(p.ASEB' .* sigma.EASE...             % ��������� � ��������� ��� ������� � ASE
        ./ (sigma.AASE + sigma.EASE));
    overlap.sum_EP  = overlap.sum_EP + sum(p.PB' .* sigma.EPB ./ ...            % ��������� � ��������� ��� �������
        (sigma.APB + sigma.EPB));
    overlap.sum_AP  = overlap.sum_AP + sum(p.PB' .* sigma.APB ./ ...            % ��������� � ����������� ��� �������
        (sigma.APB + sigma.EPB));
    overlap.sum_all = overlap.sum_all + sum(p.ASEB) + sum(p.PB);                % ����� ���� ����������� ���������
end

overlap.ES      = overlap.sum_AP ./ overlap.sum_all.^2 .* (overlap.sum_all...  % �������� ���������� (��������) ��� �������
    .* (1 - exp(-(10^(-5) ./ w_edf.S).^2)));
overlap.AS      = (overlap.sum_ES + overlap.sum_EP) ./ overlap.sum_all.^2 .*...% �������� ���������� (����������) ��� �������
    (overlap.sum_all .* (1 - exp(-(10^(-5) ./ w_edf.S).^2)));
overlap.EASE    = overlap.sum_AP ./ overlap.sum_all.^2 .* (overlap.sum_all .*...% �������� ���������� (��������) ��� ASE
    (1 - exp(-(10^(-5) ./ w_edf.ASE).^2)));
overlap.AASE    = (overlap.sum_ES + overlap.sum_EP) ./ overlap.sum_all.^2 .*... % �������� ���������� (����������) ��� ASE
    (overlap.sum_all .* (1 - exp(-(10^(-5) ./ w_edf.ASE).^2)));

g.ES            = n_sum .* sigma.ES .* overlap.ES;                              % ���. ����������� �������� (������)
g.AS            = n_sum .* sigma.AS .* overlap.AS;                              % ���. ����������� ���������� (������)
g.EASE          = n_sum .* sigma.EASE .* overlap.EASE;                          % ���. ����������� �������� (ASE)
g.AASE          = n_sum .* sigma.AASE .* overlap.AASE;                          % ���. ����������� ���������� (ASE)

alpha.PF        = (1 - exp(-(10^(-5) ./ w_edf.PF).^2)) * n_sum .* sigma.APF;    % ����������� ��������� �������� �������
alpha.PB        = (1 - exp(-(10^(-5) ./ w_edf.PB).^2)) * n_sum .* sigma.APB;    % ���������� ��������� ��������� �������

% �������� �������� ��������������

%��������� ��� ������� (dP_s/dz)
f(1:N.S,1)  = (g.ES - g.AS)' .* P(1: N.S, 1);

% ��������� ��� ��������� ASE (dP_ase/dz)
f(N.S+1:N.S+N.ASE,1)  = (g.EASE - g.AASE)' .* P(N.S+1:N.S+N.ASE, 1) + 2 * ph_const.h * ph_const.c^2 * 0.1 *...
    10^18 ./ wl.ASE'.^3 .* g.EASE';

% ��������� ��� �������� ������� (dP_p/dz)
f(N.S+N.ASE+1:N.S+N.ASE+N.PF,1)  = -alpha.PF' ./ (1 + sum(p.PF)) .* P(N.S+N.ASE+1:N.S+N.ASE+N.PF,1);

% ��������������� ��������� ��� ��������� �������
if isempty(P_in.PB) == 0
    % ��������� ��� ���������� ASE (dP_ase/dz)
    f(N.S+N.ASE+N.PF+1:N.S+N.ASE+N.PF+N.ASE,1)  = (-1) .* (g.EASE - g.AASE)' .*...
        P(N.S+N.ASE+N.PF+1: N.S+N.ASE+N.PF+N.ASE, 1) - 2 * ph_const.h * ph_const.c^2 * 0.1 * 10^18 ./ wl.ASE'.^3 ...
        .* g.EASE';
    
    % ��������� ��� ��������� ������� (dP_p/dz)
    f(N.S+N.ASE+N.PF+N.ASE+1:N.S+N.ASE+N.PF+N.ASE+N.PB,1)  = alpha.PB' ./ (1 + sum(p.PB)) .*...
        P(N.S+N.ASE+N.PF+N.ASE+1: N.S+N.ASE+N.PF+N.ASE+N.PB, 1);   
end
end