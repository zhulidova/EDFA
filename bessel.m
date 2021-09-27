%% ������ �������� ������� ������������� ������������� � �������
% ����������� �� ������� ������� ��� ������ ����� �����
function [psi, w_edf] = bessel(r_edf, Lambda, N, NA_edf, N0)
% ����������
% V_edf - ��������� ������������� ������ ��� ������� (V_edf.s), ������� (V_edf.pf � V_edf.pb) � ASE (V_edf.ase)
% U_edf � W_edf - ��������� ���������� ��� ������ ������� � ����������� �� ����� �����
% w_edf - ��������� ������� �������� ��� ������� (w_edf.s), ������� (w_edf.pf � w_edf.PB) � ASE (w_edf.ASE)

% ������ ������������� ������
V_edf.s                = 2 * pi * r_edf * NA_edf ./ Lambda.s;    % ������������� ������� �������
V_edf.pf               = 2 * pi * r_edf * NA_edf ./ Lambda.pf;   % ������������� ������� �������� �������
V_edf.ase              = 2 * pi * r_edf * NA_edf ./ Lambda.ase;  % ������������� ������� ASE

% ��������� ��� ������� ������� �������
U_edf.s                = (1 + sqrt(2)) .* V_edf.s ./ (1 + (4 + V_edf.s.^4).^0.25);
U_edf.pf               = (1 + sqrt(2)) .* V_edf.pf ./ (1 + (4 + V_edf.pf.^4).^0.25);
U_edf.ase              = (1 + sqrt(2)) .* V_edf.ase ./ (1 + (4 + V_edf.ase.^4).^0.25);
W_edf.s                = (V_edf.s.^2 - U_edf.s.^2).^0.5;
W_edf.pf               = (V_edf.pf.^2 - U_edf.pf.^2).^0.5;
W_edf.ase              = (V_edf.ase.^2 - U_edf.ase.^2).^0.5;

% ������� �������
w_edf.s                = 1.2*r_edf .* V_edf.s .* besselk(1, W_edf.s)...    % ������� ������ (������)
    .* besselj(0, U_edf.s) ./ (U_edf.s .* besselk(0, W_edf.s)); 
w_edf.pf               = 1.2*r_edf .* V_edf.pf .* besselk(1, W_edf.pf)...  % ������� ������ (�������� �������)
    .* besselj(0, U_edf.pf) ./ (U_edf.pf .* besselk(0, W_edf.pf));
w_edf.ase              = 1.2*r_edf .* V_edf.ase .* besselk(1, W_edf.ase)...% ������� ������ (ASE)
    .* besselj(0, U_edf.ase) ./ (U_edf.ase .* besselk(0, W_edf.ase));

% ������ ����������� ������� ������������� ������������� �� ���������� �� ������ �������
% ���  r < r_edf
r  = 0 : r_edf / N0 : r_edf;
n_r   = size(r,2);                           % ������ �������

psi.s(1: N.s,1: n_r)     = (besselj(0, (repmat(U_edf.s',1,n_r) .* repmat(r, N.s,1)) / r_edf)).^2; 
psi.pf(1: N.pf,1: n_r)   = (besselj(0, (repmat(U_edf.pf',1,n_r) .* repmat(r, N.pf,1)) / r_edf)).^2;
psi.ase(1: N.ase,1: n_r) = (besselj(0, (repmat(U_edf.ase',1,n_r) .* repmat(r, N.ase,1)) / r_edf)).^2;

psi.ns    = psi.s ./ pi ./ w_edf.s'.^2;      % ������������� ������� ������������� (������)
psi.npf   = psi.pf ./ pi ./ w_edf.pf'.^2;    % ������������� ������� ������������� (�������� �������)
psi.nase  = psi.ase ./ pi ./ w_edf.ase'.^2;  % ������������� ������� ������������� (ASE)

% ������ ���������� ��� ��������� �������, ���� ��� ����
if isempty(Lambda.pb) == 0
    V_edf.pb = 2 * pi * r_edf * NA_edf ./ Lambda.pb;   % ������������� ������� ��������� �������
    U_edf.pb = (1 + sqrt(2)) .* V_edf.pb ./ (1 + (4 + V_edf.pb.^4).^0.25);
    W_edf.pb = (V_edf.pb.^2 - U_edf.pb.^2).^0.5;
    w_edf.pb = r_edf .* V_edf.pb .* besselk(1, W_edf.pb)...  % ������� ������ (��������� �������)
    .* besselj(0, U_edf.pb) ./ (U_edf.pb .* besselk(0, W_edf.pb));
    psi.pb(1: N.pb,1: n_r)   = (besselj(0, (repmat(U_edf.pb',1,n_r) .* repmat(r, N.pb,1)) / r_edf)).^2;
    psi.npb            = psi.pb ./ pi ./ w_edf.pb'.^2;       % ������������� ������� ������������� (��������� �������)
end

end