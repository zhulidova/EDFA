%% ������ ������� ������������ ��������
function g = Gain(P_out, P_in)
g  = dbm(P_out)' - P_in;
end