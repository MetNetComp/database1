function [outputArg1,outputArg2] = checkResult(model)
%ANALYZERESULT この関数の概要をここに記述
%   詳細説明をここに記述
m=size(model.mets,1);

for i=1:m
%for i=4:4
    i
    list2(i,1)=0;
    list2(i,2)=0;
    
    s=sprintf('results/%d.mat',i);
    if exist(s)~=0
        load(s,'success');
        if success==1
            load(s,'gvalue');
            [GR,PR]=GRPRchecker(model,model.mets{i},gvalue);
            list2(i,1)=GR;
            list2(i,2)=PR;
            list2(i,3)=sum(cell2mat(gvalue(:,2)));
        end
    end
end

load('calculateTMGRPR.mat','list0');
list3=horzcat(list0,list2);

save('checkResult.mat');


end

