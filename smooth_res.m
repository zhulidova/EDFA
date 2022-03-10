%% Функция сглаживания спектров (для графиков)
function [x_smooth, y1_smooth, y2_smooth] = smooth_res(x,y1,y2)
x_smooth      = linspace(min(x), max(x), 25);                  % новая сетка (500 точек длин волн)
y1_smooth     = spline(x, y1, x_smooth);                        % сглаженный спектр 1              
y2_smooth     = spline(x, y2, x_smooth);                        % сглаженный спектр 2
end