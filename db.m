%% функция пересчета из раз в дБ
function d = db(x)
    d = 10 * log10(x); 
end