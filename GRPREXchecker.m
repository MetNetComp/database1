function [GR ,PR] = GRPREXchecker(model,targetMet,givenGvalue)

options=cplexoptimset('cplex');
options.mip.tolerances.integrality=10^(-12);

[model,targetRID,extype] = modelSetting(model,targetMet)

m=size(model.mets,1);
n=size(model.rxns,1);
g=size(model.genes,1);
gid=find(model.c);
pid=targetRID;
model2=model;

[grRules0] = calculateGR(model,givenGvalue);
lb2=model.lb;
ub2=model.ub;

for i=1:n
    if grRules0{i,4}==0
        lb2(i)=0;
        ub2(i)=0;
    end
end
[opt0.x, opt0.f, opt0.stat, opt0.output] = ...
    cplexlp(-model.c, [],[], model.S, zeros(m,1),lb2, ub2);
[opt0.x(gid) opt0.x(pid)]

[ex] = findExReactions(model);
s=size(ex.R,1);
k=1;kk=1;j=1;
for i=1:s
    list{i,1}=ex.R2{i};
    list{i,2}=opt0.x(ex.R(i,1));
end
[B,I]=sort(cell2mat(list(:,2)));
num_sub=size(find(B<=-0.00001),1);
num_pro=size(find(B>=0.00001),1);
for i=1:num_sub
   substrate{i,1}=list{I(i,1),1};
   substrate{i,2}=list{I(i,1),2};
end
for i=1:num_sub
   product{i,1}=list{I(s+1-i,1),1};
   product{i,2}=list{I(s+1-i,1),2};
end

file=fopen('flow.json','w');
fprintf(file,'{');
for i=1:n
   fprintf(file,'"%s": %d',model.rxns{i},opt0.x(i));
   if i~=n
       fprintf(file,', ');
   end
end
fprintf(file,'}');
fclose(file);


GR0=-opt0.f;
lb2(gid)=GR0;
ub2(gid)=GR0;
model2.c(gid)=0;
model2.c(pid)=1;
[opt1.x, opt1.f, opt1.stat, opt1.output] = ...
    cplexlp(-model2.c, [],[], model.S, zeros(m,1),lb2, ub2);
GR=opt1.x(gid) 
PR=opt1.x(pid)

save('GRPREXchecker.mat');
return;
end

