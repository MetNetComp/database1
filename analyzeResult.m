function [outputArg1,outputArg2] = analyzeResult(inputArg1,inputArg2)
%ANALYZERESULT この関数の概要をここに記述
%   詳細説明をここに記述
load('e_coli_core.mat');
model=e_coli_core;
m=size(model.mets,1);

[opt.x, opt.f, opt.stat, opt.output] = ...
    cplexlp(-model.c, [],[], model.S, zeros(m,1),model.lb, model.ub);
for i=1:m
    [model2,targetRID,extype] = modelSetting(model,model.mets{i});
    [opt.x, opt.f, opt.stat, opt.output] = ...
        cplexlp(-model2.c, [],[], model2.S, zeros(m,1),model2.lb, model2.ub);
    grid=find(model2.c);
    model2.ub(grid)=-opt.f;
    model2.lb(grid)=-opt.f;
    model2.c(grid)=0;
    model2.c(targetRID)=1;
    [opt2.x, opt2.f, opt2.stat, opt2.output] = ...
        cplexlp(-model2.c, [],[], model2.S, zeros(m,1),model2.lb, model2.ub);
   minFVA(i,1)=-opt2.f;    
end

for i=1:m
    i
    s=sprintf('results/%d.mat',i);
    load(s,'gr','pr','targetMet');
    minFVA(i,2)=gr;
    minFVA(i,3)=pr;
    minFVA(i,4)=targetMet;
end

    
save('analyzeResult.mat');


end

