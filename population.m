%% ������ ������������� �������� � ������� �������������� �������
function [n1, n2] = population(P, Lambda, sigma, psi, N, ph_const, Pin)
% ��������� 
% P_sat   - ��������� �������� ������� ��������� �������(P_sat.S), ������� (P_sat.PF � P_sat.PB) � ASE (P_sat.ASE)
% W12     - ����������� ������������ �������� � ����������� ��������� ��� �������
% W21     - ����������� ������������ �������������� �������� ��� �������
% R12     - ����������� ������������ �������� � ����������� ��������� ��� �������
% R21     - ����������� ������������ �������������� �������� ��� �������
% n1 � n2 - ������������ �������� � ������� ������

P_sat.s   = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.s .* (sigma.as + sigma.es));                  % �������� ��������� ��� �������
P_sat.ase = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.ase .* (sigma.aase + sigma.ease));            % �������� ��������� ��� ASE
P_sat.pf  = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.pf .* (sigma.apf + sigma.epf));               % �������� ��������� ��� �������� �������
% ����������� ������������ �������� � ����������� ��������� ��� ������� � ASE
W12       = sum(psi.ns .* repmat(P(1: N.s, 1),1,size(psi.ns,2)) .* repmat(sigma.as',1,size(psi.ns,2))...
    ./ ((repmat(sigma.as',1,size(psi.ns,2)) + repmat(sigma.es',1,size(psi.ns,2))) .* repmat(P_sat.s',1,size(psi.ns,2))))+...
    sum(psi.nase .* repmat(P(N.s+1: N.s+N.ase, 1),1,size(psi.nase,2)) .* repmat(sigma.aase',1,size(psi.nase,2))...
    ./ ((repmat(sigma.aase',1,size(psi.nase,2)) + repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))));
% ����������� ������������ �������������� �������� ��� ������� � ASE 
W21       = sum(psi.ns .* repmat(P(1: N.s, 1),1,size(psi.ns,2)) .* repmat(sigma.es',1,size(psi.ns,2))...
    ./ ((repmat(sigma.as',1,size(psi.ns,2)) + repmat(sigma.es',1,size(psi.ns,2))) .* repmat(P_sat.s',1,size(psi.ns,2))))+...
    sum(psi.nase .* repmat(P(N.s+1: N.s+N.ase, 1),1,size(psi.nase,2)) .* repmat(sigma.ease',1,size(psi.nase,2))...
    ./ ((repmat(sigma.aase',1,size(psi.nase,2)) + repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))));   
% ����������� ������������ �������� � ����������� ��������� ��� �������
R12    = sum(repmat(P(N.s+N.ase+1: N.s+N.ase+N.pf, 1),1,size(psi.npf,2)) .* psi.npf .* repmat(sigma.apf',1,size(psi.npf,2))...
    ./ ((repmat(sigma.apf',1,size(psi.npf,2)) + repmat(sigma.epf',1,size(psi.npf,2))) .* repmat(P_sat.pf',1,size(psi.npf,2))));
% ����������� ������������ �������������� �������� ��� �������
R21    = sum(repmat(P(N.s+N.ase+1: N.s+N.ase+N.pf, 1),1,size(psi.npf,2)) .* psi.npf .* repmat(sigma.epf',1,size(psi.npf,2))...
    ./ ((repmat(sigma.apf',1,size(psi.npf,2)) + repmat(sigma.epf',1,size(psi.npf,2))) .* repmat(P_sat.pf',1,size(psi.npf,2))));

% ��������������� ��������� ��� ��������� �������
if isempty(Pin.pb) == 0 
    P_sat.pb = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.pb .* (sigma.apb + sigma.epb));
    % ����������� ������������ �������������� �������� ��� ������� � ASE
    W21      = W21 + sum(psi.nase .* repmat(P(N.s+N.ase+N.pf+1: N.s+2*N.ase+N.pf, 1),1,size(psi.nase,2)) .*...
        repmat(sigma.ease',1,size(psi.nase,2)) ./ ((repmat(sigma.aase',1,size(psi.nase,2)) +...
        repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))));
    % ����������� ������������ �������� � ����������� ��������� ��� �������
    R12      = R12 + sum(repmat(P(N.s+2*N.ase+N.pf+1: N.s+2*N.ase+N.pf+N.pb, 1),1,size(psi.npb,2)) .*...
        psi.npb .* repmat(sigma.apb',1,size(psi.npb,2)) ./ ((repmat(sigma.apb',1,size(psi.npb,2)) +...
        repmat(sigma.epb',1,size(psi.npb,2))) .* repmat(P_sat.pb',1,size(psi.npb,2))));
    % ����������� ������������ �������������� �������� ��� �������
    R21      = R21 + sum(repmat(P(N.s+2*N.ase+N.pf+1: N.s+2*N.ase+N.pf+N.pb, 1),1,size(psi.npb,2)) .*...
        psi.npb .* repmat(sigma.epb',1,size(psi.npb,2)) ./ ((repmat(sigma.apb',1,size(psi.npb,2)) +...
        repmat(sigma.epb',1,size(psi.npb,2))) .* repmat(P_sat.pb',1,size(psi.npb,2))));
end

% ������������ �������� � ������� ������    
n2           = (W12 + R12) ./ (W12 + R12 + W21 + 1 / ph_const.tau + R21);                           % ������������ �������� ������
n1           = (1 / ph_const.tau + W21 + R21) ./ (W12 + R12 + W21 + 1 / ph_const.tau + R21);        % ������������ ������� ������

end