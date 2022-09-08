function [outputArg1,outputArg2] = trimCoreResult2()
%ANALYZERESULT この関数の概要をここに記述
%   詳細説明をここに記述
load('calculateTMGRPR.mat','list0');
load('e_coli_core.mat');
model=e_coli_core;
m=size(model.mets,1);

total=0;
success=0;
for ii=1:m
    %for i=a:b
    ii
    s=sprintf('results/trimGdel%d.mat',ii);
    if exist(s)~=0
        load(s,'size1','size2','size3','finalGRPR','time');
        if list0(ii,2)>=0.001
            %save('a.mat');
            grRatio(ii,1)=finalGRPR(1,3)/list0(ii,1);
            prRatio(ii,1)=finalGRPR(1,4)/list0(ii,2);
        else
            grRatio(ii,1)=-1;
            prRatio(ii,1)=-1;
        end
        table(ii,:)=horzcat(list0(ii,:), finalGRPR(1,3:4),grRatio(ii,1),prRatio(ii,1),size1,size2,size3);
    else
        table(ii,:)=[list0(ii,:) 0 0 0 0 0 0 0];
    end
    if table(ii,2)>=0.001
        total=total+1;
    end
    if ((table(ii,5)>=0.001) && (table(ii,6)>=0.001))
        success=success+1;
    end
end




save('trimCoreResult2.mat');


end

