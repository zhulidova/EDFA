function f = odu2(z, P, Lambda, sigma, psi, N, n_sum, ph_const, r_edf, N0)
%% �������� ���������� ��������� (����� ������)
% ����������
% r   - ������ ��������� ������� �������
% N_r - ������ ������� ��������� ������� ������� 
% n1  - ������������ ������� ������
% n2  - ������������ �������� ������
% f   - �������� ����������� � ���������� ����������

% �������
% population(P,...,r_edf) - ������ ������������� �������������� ������� � ����������� �� �������� � ����� �����

r  = 0 : r_edf / N0.r : r_edf;                                         % ������ ��������� ������� �������
N_r = size(r,2);                                                       % ������ ������� ��������� ������� �������

n  = population(P, Lambda, sigma, psi, N, ph_const,n_sum);     % ������ ������������� �������

%��������� ��� ������� (dP_s/dz)
f(1:N.s,1) = 2 * pi * n_sum * trapz(r,((repmat(sigma.es',1,N_r)-repmat(sigma.esa_eta.*sigma.as',1,N_r)) .* repmat(n.second,N.s,1)...
    - repmat(sigma.as',1,N_r) .* repmat(n.first,N.s,1) - sigma.bg).* repmat(P(1:N.s,1),1,N_r) .* psi.ns .* repmat(r,N.s,1),2);

% ��������� ��� ��������� ASE (dP_ase/dz)
f(N.s+1:N.s+N.ase,1)  = 2 * pi * n_sum * (trapz(r,((repmat(sigma.ease',1,N_r)-repmat(sigma.esa_eta.*sigma.aase',1,N_r)) .* repmat(n.second,N.ase,1)...
    .*(repmat(P(N.s+1:N.s+N.ase,1),1,N_r) + repmat(2 * ph_const.h * ph_const.c^2 * 0.1 * 10^(-9) ./...
    Lambda.ase'.^3,1,N_r)) - repmat(sigma.aase',1,N_r) .* repmat(n.first,N.ase,1) .* repmat(P(N.s+1:N.s+N.ase,1),1,N_r) - ...
    repmat(P(N.s+1:N.s+N.ase,1),1,N_r) * sigma.bg) .* psi.nase .* repmat(r,N.ase,1),2));

for i = 1: length(Lambda.pf)
    if Lambda.pf(i) > 1450* 10^(-9)
        % ��������� ��� �������� ������� (dP_p/dz) �� ����� ����� 1480 ��
        f(N.s+N.ase+i,1)  = 2 * pi * n_sum * (trapz(r,((sigma.epf(i)-sigma.apf(i)*sigma.esa_eta) .* n.second(1,:) - sigma.apf(i) .* n.first(1,:) - sigma.bg)...
            .* P(N.s+N.ase+i) .*  psi.npf(i,:) .* r(1,:)));
    elseif Lambda.pf(i) < 1100* 10^(-9)
        % ��������� ��� �������� ������� (dP_p/dz) �� ����� ����� 980 ��
        f(N.s+N.ase+i,1)  = 2 * pi * n_sum * (trapz(r,(sigma.epf(i) .* n.third(1,:) - sigma.apf(i) .* n.first(1,:) - sigma.bg) .* P(N.s+N.ase+i) .*...
            psi.npf(i,:) .* r(1,:)));
    end
end

% ��������������� ��������� ��� ��������� �������
if isempty(Lambda.pb) == 0
    % ��������� ��� ���������� ASE (dP_ase/dz)
    f(N.s+N.ase+N.pf+1:N.s+N.ase+N.pf+N.ase,1) = -2 * pi * n_sum * (trapz(r,((repmat(sigma.ease',1,N_r)-repmat(sigma.esa_eta.*sigma.aase',1,N_r)) .*...
        repmat(n.second,N.ase,1) .* (repmat(P(N.s+N.ase+N.pf+1:N.s+N.ase+N.pf+N.ase,1),1,N_r) + ...
        repmat(2 * ph_const.h * ph_const.c^2 * 0.1 * 10^(-9) ./ Lambda.ase'.^3,1,N_r)) - repmat(sigma.aase',1,N_r) .*...
        repmat(n.first,N.ase,1) .* repmat(P(N.s+N.ase+N.pf+1:N.s+N.ase+N.pf+N.ase,1),1,N_r) - repmat(P(N.s+N.ase+N.pf+1:N.s+N.ase+N.pf+N.ase,1),1,N_r)*sigma.bg)...
        .* psi.nase .* repmat(r,N.ase,1),2));
    for i = 1: length(Lambda.pb)
        if Lambda.pb(i) > 1450* 10^(-9)
            % ��������� ��� ��������� ������� (dP_p/dz) �� ����� ����� 1480 ��
            f(N.s+N.ase+N.pf+N.ase+i,1)  = -2 * pi * n_sum * (trapz(r,((sigma.epb(i)-sigma.apb(i)*sigma.esa_eta) .* n.second(1,:) - sigma.apb(i) .* n.first(1,:) - sigma.bg)...
                .* P(N.s+N.ase+N.pf+N.ase+i) .*  psi.npb(i,:) .* r(1,:)));
        elseif Lambda.pb(i) < 1100* 10^(-9)
            % ��������� ��� ��������� ������� (dP_p/dz) �� ����� ����� 980 ��
            f(N.s+N.ase+N.pf+N.ase+i,1)  = -2 * pi * n_sum * (trapz(r,(sigma.epb(i) .* n.third(1,:) - sigma.apb(i) .* n.first(1,:) - sigma.bg) .* P(N.s+N.ase+N.pf+N.ase+i) .*...
                psi.npb(i,:) .* r(1,:)));
        end
    end
end
end