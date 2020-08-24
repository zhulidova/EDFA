%% Расчет спектра коэффициента усиления
function Gain = gain(P_out, P_in)
Gain  = dbm(P_out) - dbm(P_in);
end