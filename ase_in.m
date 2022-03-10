%% ������� ������� ������� �������� ASE
%(Desurvire, appendix R)
% P_sat_P   - ������ �������� �������� ��������� ��� �������
% P_sat_ASE - ������ �������� �������� ��������� ��� ������� ASE
% n_sp      - ������ ������������� ���������� �������������, ����� ������ ������� ASE
% p_ase     - �������� ASE � ��
% P_ase     - �������� ASE � ���

function P_ase = ase_in(Lambda, ph_const, sigma, Pin, n_sum, psi, N, r_edf, N0, L)

r  = 0 : r_edf / N0.r : r_edf;  % ������ ��������� ������� �������
N_r = size(r,2);

% ���� ��� ��� ��������� ASE ��� ��������� �����������, ��������� �������
% ���������������� ��� �������� � ������������� ����������� �������
if isempty(Pin.asef) == 0
    Lambda.pf   = Lambda.pb;
    Pin.pf      = Pin.pb;
    sigma.apf   = sigma.apb;
    sigma.epf   = sigma.apb;
    psi.pf      = psi.pb;
    psi.npf     = psi.npb;
    N.pf        = N.pb;
    sigma.bgpf  = sigma.bgpb;
end

% ��� ������� ASE ���� ������ ����� �������, ��������
Lambda.pb   = [];
Pin.pb      = [];
sigma.apb   = [];
sigma.epb   = [];
psi.pb      = [];
psi.npb     = [];
N.pb        = [];
sigma.bgpb  = [];
Pin.ases    = zeros(1,N.s);
Pin.asef    = zeros(1,N.ase);

% ������ �������� ASE � ����� �����
[z, Pout] = ode45(@(z,P) odu2(z, P, Lambda, sigma, psi, N, n_sum, ph_const, r_edf, N0), [0: L/N0.z: L],[Pin.s Pin.asef Pin.pf]);

P.ase     = Pout(:, 1 + N.s: N.ase + N.s)';  % ������������� ��������� ��������

for i = 1:size(P.ase,2)
    Gain_z(:,i)  = P.ase(:,i) - Pin.asef';
    n = population(Pout(i, :)', Lambda, sigma, psi, N, ph_const,n_sum);
    a(:,i) = 2 * pi * n_sum * trapz(r,(repmat(sigma.ease',1,N_r) .* repmat(n.second,N.ase,1)).* repmat(P.ase(:,i),1,N_r) .* psi.nase .* repmat(r,N.ase,1),2);
    d_n_eq(:,i) = a(:,i)./Gain_z(:,i);
end

n_eq  = a(:,size(a,2)) ./ Gain_z(:,size(Gain_z,2)); % ������ ������������ ������������� ����

P_ase = n_eq' .* ph_const.c * 0.1 * 10^(-9) ./ Lambda.ase.^3 .* ph_const.h * ph_const.c; % �������� ASE � ��

end