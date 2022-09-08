function [] = test()
initCobraToolbox
load('e_coli_core.mat');
model=e_coli_core;
[vg gr pr it success]=gDel_minRN(model,'succ_e',10,0.001,0.001);

save('test.mat');
end

