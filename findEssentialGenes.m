function [GR ,PR] = findEssentialGenes(model)

options=cplexoptimset('cplex');
options.mip.tolerances.integrality=10^(-12);

m=size(model.mets,1);
n=size(model.rxns,1);
g=size(model.genes,1);
gid=find(model.c);

for j=1:g
    gvalue0{j,1}=model.genes{j};
    gvalue0{j,2}=1;
end

for i=1:g
    i
    model2=model;
    gvalue=gvalue0;
    gvalue{i,2}=0;
    
    [grRules0] = calculateGR(model,gvalue);
    
    for j=1:n
        if grRules0{j,4}==0
            model2.lb(j)=0;
            model2.ub(j)=0;
        end
    end
    [opt0.x, opt0.f, opt0.stat, opt0.output] = ...
        cplexlp(-model.c,[],[],model.S,zeros(m,1),model2.lb,model2.ub);
    if opt0.stat>0
        if opt0.x(gid)>=0.001
            essential(i,1)=0;
        else
            essential(i,1)=1;
        end
    else
        essential(i,1)=1;
    end
    
end
save('findEssentialGenes.mat');
return;
end

