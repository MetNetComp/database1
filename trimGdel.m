function [gvalue, finalGRPR, size1, size2, size3,time] = trimGdel(model,targetMet,givenGvalue,i)
tic;
options=cplexoptimset('cplex');
options.mip.tolerances.integrality=10^(-12);
sss=sprintf('results/trimGdel%d.mat',i);
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

GR0=-opt0.f;
lb2(gid)=GR0;
ub2(gid)=GR0;
model2.c(gid)=0;
model2.c(pid)=1;
[opt1.x, opt1.f, opt1.stat, opt1.output] = ...
    cplexlp(-model2.c, [],[], model.S, zeros(m,1),lb2, ub2);
[opt1.x(gid) opt1.x(pid)]
GRLB=opt1.x(gid);
PRLB=opt1.x(pid);
[term,ng,nt,nr,nko,reactionKO,reactionKO2term] = readGeneRules(model);
[f,intcon,A,b,Aeq,beq,lb,ub,xname] = geneReactionMILP(model,term,ng,nt,nr,nko,reactionKO);

lp.Aeq=Aeq;
lp.beq=[zeros(size(lp.Aeq,1),1)];
j=1;
for i=1:size(model.grRules,1)
    if isempty(model.grRules{i,:})==0
        ind(1,j)=i;
        j=j+1;
    end
end
z1=-diag(model.ub);
z2=diag(model.lb);
z3=eye(n);
lp.A=A;
lp.b=b;
lp.lb=lb;
lp.ub=ub;
for i=1:ng
    if givenGvalue{i,2}==1
        lp.lb(i,1)=1;
        lp.ub(i,1)=1;
    end
end

[grRules0] = calculateGR(model,givenGvalue);
j=1;
for i=1:n
    if isempty(model.grRules{i,1})==0
        lp.lb(ng+nt+j,1)=grRules0{i,4};
        lp.ub(ng+nt+j,1)=grRules0{i,4};
        j=j+1;
    end
end
lp.f=[-ones(ng,1); zeros(nt,1); zeros(nko,1)];
for i=1:n
    s2=repelem('B',ng+nt+nko);
    lp.ctype=sprintf('%s%s',s2);
end

[opt.x, opt.f, opt.stat, opt.output] = ...
    cplexmilp(lp.f, lp.A, lp.b, lp.Aeq, lp.beq,[],[],[],lp.lb, lp.ub,lp.ctype,[],options);

gvalue=givenGvalue;
if opt.stat>0
    for i=1:ng
        vg(i,1)=opt.x(i);
        gvalue{i,2}=opt.x(i);
    end
    for i=1:nt
        vt(i,1)=opt.x(ng+i);
    end
    for i=1:nko
        vko(i,1)=opt.x(ng+nt+i);
    end
end
gvalue0=gvalue;
trimmed=1;
model2=model;
grprlist(1,1)=opt1.x(gid);
grprlist(1,2)=opt1.x(pid);
grprlist(1,3)=opt1.x(gid);
grprlist(1,4)=opt1.x(pid);
k=2;

while trimmed==1;
    trimmed=0;
    for i=1:ng
        %i
        if gvalue{i,2}==0
            gvalue{i,2}=1;
            [grRules2] = calculateGR(model,gvalue);
            lb2=model.lb;
            ub2=model.ub;

            for j=1:n
                if grRules2{j,4}==0
                    lb2(j)=0;
                    ub2(j)=0;
                end
            end
            [opt2.x, opt2.f, opt2.stat, opt2.output] = ...
                cplexlp(-model.c, [],[], model.S, zeros(m,1),lb2, ub2);
            grprlist(k,1)=opt2.x(gid);
            grprlist(k,2)=opt2.x(pid);
            GR2=-opt2.f;
            lb2(gid)=GR2;
            ub2(gid)=GR2;
            model2.c(gid)=0;
            model2.c(pid)=1;
            [opt3.x, opt3.f, opt3.stat, opt3.output] = ...
                cplexlp(model2.c, [],[], model.S, zeros(m,1),lb2, ub2);
            grprlist(k,3)=opt3.x(gid);
            grprlist(k,4)=opt3.x(pid);
            if  ((opt3.x(gid) < 0.999 * GRLB) || (opt3.x(pid) < 0.999*PRLB))
                gvalue{i,2}=0;
                grprlist(k,:)=grprlist(k-1,:);
                %s=sprintf('%d is not trimmed.',i);display(s)
            else
                trimmed=1;
                %s=sprintf('%d is trimmed.',i);display(s)
            end
            finalGRPR=[GRLB PRLB grprlist(k,3) grprlist(k,4)];
            k=k+1;
        end
    end
end

gvalueList=horzcat(givenGvalue(:,2),gvalue0(:,2),gvalue(:,2));
size1=size(find(cell2mat(givenGvalue(:,2))==0),1);
size2=size(find(cell2mat(gvalue0(:,2))==0),1);
size3=size(find(cell2mat(gvalue(:,2))==0),1);

finalGRPR;
time=toc;
save(sss);
return;
end

