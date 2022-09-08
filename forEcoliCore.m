function [outputArg1,outputArg2] = forEcoliCore()
%example2 calculates the gene deletion strategy for growth coupling
%for biotin in iML1515.
%
% Apr. 23, 2021  Takeyuki TAMURA
%
load('e_coli_core.mat');
model=e_coli_core;
m=size(model.mets,1);
%for i=46:46
for i=1:m
    model=e_coli_core;
    [gr pr it success]=gDel_minRN(model,model.mets{i},i,10,0.001,0.001);
end

save('forEcoliCore.mat');
end

