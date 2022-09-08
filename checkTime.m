function [outputArg1,outputArg2] = checkTime(model)
%ANALYZERESULT この関数の概要をここに記述
%   詳細説明をここに記述
m=size(model.mets,1);

for i=1:m
%for i=4:4
    i
    list(i,1)=0;
 
    s=sprintf('results/%d.mat',i);
    if exist(s)~=0
        load(s,'success','time');
        if success==1
            list(i,1)=time;
        end
    end
end



save('checkTime.mat');


end

