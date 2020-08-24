function f = odu2_unsaturated_gain_regime(z, P, Lambda, sigma, N,w_edf, n_sum, ph_const, P_in)
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

P_sat.s             = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.s.^2)...   % �������� ��������� (������)
    ./ (Lambda.s .* (sigma.as + sigma.es) * ph_const.tau);
P_sat.pf            = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.pf.^2)...  % �������� ��������� (�������� �������)
    ./ (Lambda.pf .* (sigma.apf + sigma.epf) * ph_const.tau);                                   
P_sat.pb            = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.pb.^2)...  % �������� ��������� (��������� �������)
    ./ (Lambda.pb .* (sigma.apb + sigma.epb) * ph_const.tau);                                   
P_sat.ase           = (ph_const.h * ph_const.c * 10^25 * pi .* w_edf.ase.^2)... % �������� ��������� (ASE)
    ./ (Lambda.ase .* (sigma.aase + sigma.ease) * ph_const.tau);

p.s                 = P(1: N.s, 1) ./ P_sat.s';                                 % ����������� �������� �������
p.asef              = P(N.s+1:N.s+N.ase, 1) ./ P_sat.ase';                      % ����������� �������� ��������� ASE
p.pf                = P(N.s+N.ase+1:N.s+N.ase+N.pf,1) ./ P_sat.pf';             % ����������� �������� �������� �������
overlap.sum_es      = sum(p.s' .* sigma.es ./ (sigma.as + sigma.es)) +...       % ��������� � ��������� ��� ������� � ASE
    sum(p.asef' .* sigma.ease ./ (sigma.aase + sigma.ease));
overlap.sum_all     = sum(p.s) + sum(p.asef) + sum(p.pf);                       % ����� ���� ����������� ���������
overlap.sum_ep      = sum(p.pf' .* sigma.epf ./ (sigma.apf + sigma.epf));       % ��������� � ��������� ��� �������
overlap.sum_ap      = sum(p.pf' .* sigma.apf ./ (sigma.apf + sigma.epf));       % ��������� � ����������� ��� �������

if isempty(P_in.pb) == 0
    p.aseb          = P(N.s+N.ase+N.pf+1:N.s+2*N.ase+N.pf) ./ P_sat.ase';   % ����������� �������� ��������� ASE
    p.pb            = P(N.s+2*N.ase+N.pf+1:N.s+2*N.ase+N.pf+N.pb) ...   % ����������� �������� ��������� �������
        ./ P_sat.pb';
    overlap.sum_es  = overlap.sum_es + sum(p.aseb' .* sigma.ease...             % ��������� � ��������� ��� ������� � ASE
        ./ (sigma.aase + sigma.ease));
    overlap.sum_ep  = overlap.sum_ep + sum(p.pb' .* sigma.epb ./ ...            % ��������� � ��������� ��� �������
        (sigma.apb + sigma.epb));
    overlap.sum_ap  = overlap.sum_ap + sum(p.pb' .* sigma.apb ./ ...            % ��������� � ����������� ��� �������
        (sigma.apb + sigma.epb));
    overlap.sum_all = overlap.sum_all + sum(p.aseb) + sum(p.pb);                % ����� ���� ����������� ���������
end

overlap.es      = overlap.sum_ap ./ overlap.sum_all.^2 .* (overlap.sum_all...   % �������� ���������� (��������) ��� �������
    .* (1 - exp(-(10^(-5) ./ w_edf.s).^2)));
overlap.as      = (overlap.sum_es + overlap.sum_ep) ./ overlap.sum_all.^2 .*... % �������� ���������� (����������) ��� �������
    (overlap.sum_all .* (1 - exp(-(10^(-5) ./ w_edf.s).^2)));
overlap.ease    = overlap.sum_ap ./ overlap.sum_all.^2 .* (overlap.sum_all .*...% �������� ���������� (��������) ��� ASE
    (1 - exp(-(10^(-5) ./ w_edf.ase).^2)));
overlap.aase    = (overlap.sum_es + overlap.sum_ep) ./ overlap.sum_all.^2 .*... % �������� ���������� (����������) ��� ASE
    (overlap.sum_all .* (1 - exp(-(10^(-5) ./ w_edf.ase).^2)));

g.es            = n_sum .* sigma.es .* overlap.es;                              % ���. ����������� �������� (������)
g.as            = n_sum .* sigma.as .* overlap.as;                              % ���. ����������� ���������� (������)
g.ease          = n_sum .* sigma.ease .* overlap.ease;                          % ���. ����������� �������� (ASE)
g.aase          = n_sum .* sigma.aase .* overlap.aase;                          % ���. ����������� ���������� (ASE)

alpha.pf        = (1 - exp(-(10^(-5) ./ w_edf.pf).^2)) * n_sum .* sigma.apf;    % ����������� ��������� �������� �������
alpha.pb        = (1 - exp(-(10^(-5) ./ w_edf.pb).^2)) * n_sum .* sigma.apb;    % ���������� ��������� ��������� �������

% �������� �������� ��������������

%��������� ��� ������� (dP_s/dz)
f(1:N.s,1)  = (g.es - g.as)' .* P(1: N.s, 1);

% ��������� ��� ��������� ASE (dP_ase/dz)
f(N.s+1:N.s+N.ase,1)  = (g.ease - g.aase)' .* P(N.s+1:N.s+N.ase, 1) + 2 * ph_const.h * ph_const.c^2 * 0.1 *...
    10^(-9) ./ Lambda.ase'.^3 .* g.ease';

% ��������� ��� �������� ������� (dP_p/dz)
f(N.s+N.ase+1:N.s+N.ase+N.pf,1)  = -alpha.pf' ./ (1 + sum(p.pf)) .* P(N.s+N.ase+1:N.s+N.ase+N.pf,1);

% ��������������� ��������� ��� ��������� �������
if isempty(P_in.pb) == 0
    % ��������� ��� ���������� ASE (dP_ase/dz)
    f(N.s+N.ase+N.pf+1:N.s+2*N.ase+N.pf,1)  = (-1) .* (g.ease - g.aase)' .*...
        P(N.s+N.ase+N.pf+1: N.s+2*N.ase+N.pf, 1) - 2 * ph_const.h * ph_const.c^2 * 0.1 * 10^(-9) ./ Lambda.ase'.^3 ...
        .* g.ease';
    
    % ��������� ��� ��������� ������� (dP_p/dz)
    f(N.s+2*N.ase+N.pf+1:N.s+2*N.ase+N.pf+N.pb,1)  = alpha.pb' ./ (1 + sum(p.pb)) .*...
        P(N.s+2*N.ase+N.pf+1: N.s+2*N.ase+N.pf+N.pb, 1);   
end
end