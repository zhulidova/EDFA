%% EFDA slow
% Array of Lambda(i) for diskretisation ASE
function Lambda_array=ase_diskr(Nase, Lambda_range)
% Lambda_array=[Lambda_range(1)+(Lambda_range(1)-Lambda_range(2))/2/Nase];
% for i=1:(Nase-1)
%     Lambda_array=[Lambda_array Lambda_range(1)+(Lambda_range(2)-Lambda_range(1))/Nase*(i+0.5)];
% end
     del = (Lambda_range(2)-Lambda_range(1))/(Nase-1);
     Lambda_array=Lambda_range(1):del:Lambda_range(2);
end

