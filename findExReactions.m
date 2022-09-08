function [ex] = findExReactions(model)
%FINDEXREACTIONS この関数の概要をここに記述
%   詳細説明をここに記述
m=size(model.mets,1);
n=size(model.rxns,1);

k=1;
for i=1:n
    if size(find(model.S(:,i)),1)==1
       ex.R(k,1)=i;
       ex.R2{k,1}=model.rxns{i};
       ex.met(k,1)=find(model.S(:,i));
       ex.met2{k,1}=model.mets{ex.met(k,1)};
       %c1=count(model.metFormulas{ex.met(k,1)},'N');
       %c2=count(model.metFormulas{ex.met(k,1)},'Ni');
       %c3=count(model.metFormulas{ex.met(k,1)},'Ne');
       %c4=count(model.metFormulas{ex.met(k,1)},'Na');
       %ex.nitro(k,1)=c1-c2-c3-c4;
       k=k+1;
    end
end

save('findExReactions.mat');
end

