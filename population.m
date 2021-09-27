%% ������ ������������� �������� � ������� �������������� �������
function n = population(P, Lambda, sigma, psi, N, ph_const,n_sum)
% ��������� 
% P_sat   - ��������� �������� ������� ��������� �������(P_sat.s), ������� (P_sat.pf � P_sat.pb) � ASE (P_sat.ase)
% W12     - �������� ������������ �������� � ����������� ��������� ��� �������
% W21     - �������� ������������ �������������� �������� ��� �������
% R12     - �������� ������������ �������� � ����������� ��������� ��� �������
% R21     - �������� ������������ �������������� �������� ��� �������
% R24     - �������� ������������ �������� � ����������� � ������ I(13/2)
% �� I(9/2), ����� ����� 1680 � ������������� � (11/2) �� F(7/2), ����� ����� 980
% n1      - ������ ������������� ������ 1 : I(15/2), ��� ������ ��������� � ����������
% n2      - ������ ������������� ������ 2 : I(13/2), ��� ������ ��������� � ����������
% n3      - ������ ������������� ������ 3 : I(11/2), ��� ������ ��������� � ����������
% n4      - ������ ������������� ������ 4 : I(19/2), ��� ������ ��������� � ����������

P_sat.s   = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.s .* (sigma.as + sigma.es));       % �������� ��������� ��� �������
P_sat.ase = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.ase .* (sigma.aase + sigma.ease)); % �������� ��������� ��� ASE
P_sat.pf  = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.pf .* (sigma.apf + sigma.epf));    % �������� ��������� ��� �������� �������

R12 = 0;
R21 = 0;
R13 = 0;
R31 = 0;
R24 = 0;
n.first  = zeros(1,size(psi.ns,2));
n.second = zeros(1,size(psi.ns,2));
n.third  = zeros(1,size(psi.ns,2));
n.fourth = zeros(1,size(psi.ns,2));

%% ������� ��������� ���������
% �������� ������������ �������� � ����������� ��������� ��� ������� � ASE
W12       = sum(psi.ns .* repmat(P(1: N.s, 1),1,size(psi.ns,2)) .* repmat(sigma.as',1,size(psi.ns,2))...
    ./ ((repmat(sigma.as',1,size(psi.ns,2)) + repmat(sigma.es',1,size(psi.ns,2))) .* repmat(P_sat.s',1,size(psi.ns,2))))+...
    sum(psi.nase .* repmat(P(N.s+1: N.s+N.ase, 1),1,size(psi.nase,2)) .* repmat(sigma.aase',1,size(psi.nase,2))...
    ./ ((repmat(sigma.aase',1,size(psi.nase,2)) + repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))),1);

% �������� ������������ �������������� �������� ��� ������� � ASE 
W21       = sum(psi.ns .* repmat(P(1: N.s, 1),1,size(psi.ns,2)) .* repmat(sigma.es',1,size(psi.ns,2))...
    ./ ((repmat(sigma.as',1,size(psi.ns,2)) + repmat(sigma.es',1,size(psi.ns,2))) .* repmat(P_sat.s',1,size(psi.ns,2))))+...
    sum(psi.nase .* repmat(P(N.s+1: N.s+N.ase, 1),1,size(psi.nase,2)) .* repmat(sigma.ease',1,size(psi.nase,2))...
    ./ ((repmat(sigma.aase',1,size(psi.nase,2)) + repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))),1); 

R24  = sum(psi.ns .* repmat(P(1: N.s, 1),1,size(psi.ns,2)) .* repmat(sigma.esa_eta.*sigma.as',1,size(psi.ns,2)) ./...
    (ph_const.h * ph_const.c * 10^25).* repmat(Lambda.s',1,size(psi.ns,2)),1);

R24  = R24 + sum(psi.nase .* repmat(P(N.s+1: N.s+N.ase, 1),1,size(psi.nase,2)) .* repmat(sigma.esa_eta.*sigma.aase',1,size(psi.nase,2))...
    ./ (ph_const.h * ph_const.c * 10^25).*repmat(Lambda.ase',1,size(psi.nase,2)),1);


for i = 1:length(Lambda.pf)
if Lambda.pf(i) > 1450* 10^(-9) % ��� ������� � ������ ����� 1480 ��
    
    % �������� ������������ �������� � ����������� ��������� ��� �������
    R12  = R12 + psi.npf(i,:) .* P(N.s+N.ase+i, 1) .* sigma.apf(i) ./ (sigma.apf(i) + sigma.epf(i)) ./ P_sat.pf(i);
    
    % �������� ������������ �������������� �������� ��� �������
    R21  = R21 + psi.npf(i,:) .* P(N.s+N.ase+i, 1) .* sigma.epf(i) ./ (sigma.apf(i) + sigma.epf(i)) ./ P_sat.pf(i);
    
    % ���� ������� ESA ��� ������� �� ����� ����� 1480 �� (��� 980 �� ������ ���� ������� �������� �������)
    R24  = R24 + psi.npf(i,:) .* P(N.s+N.ase+i, 1) .* sigma.esa_eta*sigma.apf(i) ./ (ph_const.h * ph_const.c * 10^25).*Lambda.pf(i);


else                            % ��� ������� � ������ ����� 980 ��
    
    % �������� ������������ �������� � ����������� ��������� ��� �������
    R13  = R13 + psi.npf(i,:) .* P(N.s+N.ase+i, 1) .* sigma.apf(i) ./ (sigma.apf(i) + sigma.epf(i)) ./ P_sat.pf(i);
    
    % �������� ������������ �������������� �������� ��� �������
    R31  = R31 + psi.npf(i,:) .* P(N.s+N.ase+i, 1) .* sigma.epf(i) ./ (sigma.apf(i) + sigma.epf(i)) ./ P_sat.pf(i);
end
end
% ��������������� ��������� ��� ��������� �������
if isempty(Lambda.pb) == 0
    P_sat.pb = (ph_const.h * ph_const.c * 10^25) ./ (Lambda.pb .* (sigma.apb + sigma.epb));
    
    % �������� ������������ �������������� �������� ��� ������� � ASE
    W21   = W21 + sum(psi.nase .* repmat(P(N.s+N.ase+N.pf+1: N.s+2*N.ase+N.pf, 1),1,size(psi.nase,2)) .*...
        repmat(sigma.ease',1,size(psi.nase,2)) ./ ((repmat(sigma.aase',1,size(psi.nase,2)) +...
        repmat(sigma.ease',1,size(psi.nase,2))) .* repmat(P_sat.ase',1,size(psi.nase,2))),1);
    
    R24  = R24 + sum(psi.nase .* repmat(P(N.s+N.ase+N.pf+1: N.s+2*N.ase+N.pf, 1),1,size(psi.nase,2)) .*...
        repmat(sigma.esa_eta.*sigma.aase',1,size(psi.nase,2)) ./ (ph_const.h * ph_const.c * 10^25).* ...
        repmat(Lambda.ase',1,size(psi.nase,2)),1);

    for i = 1:length(Lambda.pb)
        if Lambda.pb(i) > 1450* 10^(-9)  % ��� ������� � ������ ����� 1480 ��          
            
            % �������� ������������ �������� � ����������� ��������� ��� �������
            R12  = R12 + psi.npb(i,:) .* P(N.s+2*N.ase+N.pf+i,1) .*  sigma.apb(i) ./ (sigma.apb(i) + sigma.epb(i)) ./ P_sat.pb(i);
            
            % �������� ������������ �������������� �������� ��� �������
            R21  = R21 + psi.npb(i,:) .* P(N.s+2*N.ase+N.pf+i,1) .*  sigma.epb(i) ./ (sigma.apb(i) + sigma.epb(i)) ./ P_sat.pb(i);
            
            % ���� ������� ESA ��� ������� �� ����� ����� 1480 �� (��� 980 �� ������ ���� ������� �������� �������)
            R24  = R24 + psi.npb(i,:) .* P(N.s+2*N.ase+N.pf+i,1) .* sigma.esa_eta*sigma.apb(i) ./ (ph_const.h * ph_const.c * 10^25) .* Lambda.pb(i);
        else                             % ��� ������� � ������ ����� 980 ��
            
            % �������� ������������ �������� � ����������� ��������� ��� �������
            R13   = R13 + psi.npb(i,:) .* P(N.s+2*N.ase+N.pf+i) .*  sigma.apb(i) ./ (sigma.apb(i) + sigma.epb(i)) ./ P_sat.pb(i);
            
            % �������� ������������ �������������� �������� ��� �������
            R31   = R31 + psi.npb(i,:) .* P(N.s+2*N.ase+N.pf+i) .*  sigma.epb(i) ./ (sigma.apb(i) + sigma.epb(i)) ./ P_sat.pb(i);
        end
    end
end

if isempty(Lambda.pb) == 1
    Lambda.pb = Lambda.pf;
end
%% ������ ������������� ���������������� ������� (���������� ������������� ������� ����� ����� ����������� � ESA)

% ���� ������������� �����������
C24 = (2.65 * n_sum * 0.1 + 3.38) * 10; % ���� ������������� ����������� � ������ 2 �� 4 (1-�� ���������:
                                        % ����������� ������������� �����������, ������� �� ������������;
                                        % 2-�� ���������: ����������� ��������� ����� ��������� I(11/2)-I(15/2) (980 ��)
                                        % � ��������� I(11/2)-I(13/2) (���������������� �������))

% ��� ������� � ������ ����� 1480 �� 
if min(Lambda.pf) > 1450* 10^(-9) || min(Lambda.pb) > 1450* 10^(-9)
    
    % ������������ ��� ������������ n2
    D = (1 + W21./(W12+R12) + 1./(ph_const.tau2.*(W12+R12)) + R21./(W12+R12) + (ph_const.tau3+ph_const.tau4)*R24).^2 +...
    4.*(C24./(W12+R12) + (ph_const.tau3+ph_const.tau4)*C24);

    % ������ ������������� ������ 2 : I(13/2), ��� ������ ��������� � ����������
    n.second = ((-1).*(1 + W21./(W12+R12) + 1./(ph_const.tau2.*(W12+R12)) + R21./(W12+R12) + (ph_const.tau3+ph_const.tau4)*R24) + sqrt(D))...
    ./ (2*(C24./(W12+R12) + (ph_const.tau3+ph_const.tau4)*C24));

    % ������ ������������� ������ 1 : I(15/2), ��� ������ ��������� � ����������
    n.first = (n.second.*W21 + n.second/ph_const.tau2 + C24.*n.second.^2 + R21.*n.second) ./ (W12+R12);

    % ������ ������������� ������ 4 : I(19/2), ��� ������ ��������� � ����������
    n.fourth = ph_const.tau4 * (C24*n.second.^2 + n.second.*R24);

    % ������ ������������� ������ 3 : I(11/2), ��� ������ ��������� � ����������
    n.third = ph_const.tau3 / ph_const.tau4 .* n.fourth;

% ��� ������� � ������ ����� 980 ��    
elseif max(Lambda.pf) < 1100* 10^(-9) || max(Lambda.pb) < 1100* 10^(-9) 
    
    % ������������ ��� ������������ n2
    D = (1/ph_const.tau2./(R13+W12) + W21./(R13+W12) + (R31./(R13+W12) + 1) * ph_const.tau3/ph_const.tau2 .* (R13 + W21*ph_const.tau2.*R13)...
        ./ (R13+W12+R31.*W12*ph_const.tau3) + 1).^2 + 4 * (C24./(R13+W12) + (R31./(R13+W12) + 1) * ph_const.tau3 * C24.*(2.*R13+W12)...
        ./ (R13+W12+R31.*W12*ph_const.tau3) + ph_const.tau4 * C24);
    
    % ������ ������������� ������ 2 : I(13/2), ��� ������ ��������� � ����������
    n.second = ((-1)*(1/ph_const.tau2./(R13+W12) + W21./(R13+W12) + (R31./(R13+W12) + 1)*ph_const.tau3/ph_const.tau2 .* (R13+W21*ph_const.tau2.*R13)...
        ./ (R13+W12+R31.*W12*ph_const.tau3) + 1) + sqrt(D))./2./(C24./(R13+W12) + (R31./(R13+W12) + 1)*ph_const.tau3 * C24.*(2.*R13+W12)...
        ./ (R13+W12+R31.*W12*ph_const.tau3) + ph_const.tau4 * C24);
    
    % ������ ������������� ������ 4 : I(19/2), ��� ������ ��������� � ����������
    n.fourth = C24 .* n.second.^2 * ph_const.tau4;
    
    % ������ ������������� ������ 3 : I(11/2), ��� ������ ��������� � ����������
    n.third  = ph_const.tau3/ph_const.tau2 .* (C24 .* n.second*ph_const.tau2.*(2.*R13+W12) + n.second.*R13 + W21.*n.second*ph_const.tau2.*R13)...
        ./ (R13+W12+R31.*W12*ph_const.tau3);
    
    % ������ ������������� ������ 1 : I(15/2), ��� ������ ��������� � ����������
    n.first = (n.third/ph_const.tau3 - n.fourth/ph_const.tau4 + R31.*n.third) ./ R13;

else
    % ������������ ��� ������������ n2
    D = (-R24./R13 + (1./R13/ph_const.tau3+R31./R13+1) * ph_const.tau3 .* ((R13+W12+R12).*R24 + R13/ph_const.tau2 + R13.*W21 + R13.*R21) ./...
        (R13+W12+R12+R31.*W12*ph_const.tau3+R31.*R12*ph_const.tau3) + 1 + ph_const.tau4.*R24).^2 + 4.*(-C24./R13+(1./R13/ph_const.tau3+R31./R13+1) *...
        ph_const.tau3 .* ((R13+W12+R12)*C24+R13*C24) ./ (R13+W12+R12+R31.*W12*ph_const.tau3+R31.*R12*ph_const.tau3) + ph_const.tau4*C24);
    
    % ������ ������������� ������ 2 : I(13/2), ��� ������ ��������� � ����������
    n.second = (-(-R24./R13+(1./R13/ph_const.tau3+R31./R13+1) * ph_const.tau3 .* ((R13+W12+R12).*R24+R13/ph_const.tau2+R13.*W21+R13.*R21) ./...
        (R13+W12+R12+R31.*W12*ph_const.tau3+R31.*R12*ph_const.tau3) + 1 + ph_const.tau4.*R24) + sqrt(D)) / 2. /(-C24./R13+(1./R13/ph_const.tau3+R31./R13+1) *...
        ph_const.tau3.*((R13+W12+R12)*C24+R13*C24) ./ (R13+W12+R12+R31.*W12*ph_const.tau3+R31.*R12*ph_const.tau3) + ph_const.tau4*C24);
    
    % ������ ������������� ������ 4 : I(19/2), ��� ������ ��������� � ����������
    n.fourth = ph_const.tau4 .* (C24 .* n.second.^2 + R24.*n.second);
    
    % ������ ������������� ������ 3 : I(11/2), ��� ������ ��������� � ����������
    n.third  = ph_const.tau3.*((C24.*n.second.^2+R24.*n.second).*(R13+W12+R12)+R13.*(n.second/ph_const.tau2+C24.*n.second.^2+W21.*n.second+R21.*n.second))...
        ./(R13+W12+R12+R31.*W12*ph_const.tau3+R31.*R12*ph_const.tau3);
    
    % ������ ������������� ������ 1 : I(15/2), ��� ������ ��������� � ����������
    n.first = (R21.*n.second + n.second/ph_const.tau2 - n.third/ph_const.tau3 + 2*C24.*n.second.^2 + W21.*n.second + R24 .*n.second)./(R12+W12);
end

% n2           = (W12 + R12) ./ (W12 + R12 + W21 + 1 / ph_const.tau + R21);                           % ������������ �������� ������
% n1           = (1 / ph_const.tau + W21 + R21) ./ (W12 + R12 + W21 + 1 / ph_const.tau + R21);        % ������������ ������� ������
end