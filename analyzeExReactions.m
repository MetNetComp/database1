function [outputArg1,outputArg2] = analyzeExReactions(model)
%ANALYZEEXREACTIONS この関数の概要をここに記述
%   詳細説明をここに記述
opt0 = optimizeCbModel(model);
[ex] = findExReactions(model);
s=size(ex.R,1);
k=1;kk=1;j=1;
for i=1:s
    list{i,1}=ex.R2{i};
    list{i,2}=opt0.x(ex.R(i,1));
end


save('analyzeExReactions.mat');
end

