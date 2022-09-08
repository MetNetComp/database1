function [outputArg1,outputArg2] = makeEscherInfo(model,targetMet,gvalue)
%MAKEESCHERINFO この関数の概要をここに記述
%   詳細説明をここに記述
model2=modelSetting(model,targetMet);
n=size(model2.rxns,1);
[grRules] = calculateGR(model2,gvalue);
 
for i=1:n
    if grRules{i,4}==0
       model2.ub(i)=0;
       model2.lb(i)=0;
    end
end
opt2=optimizeCbModel(model2);

file=fopen('flow.json','w');
fprintf(file,'{');
for i=1:n

   fprintf(file,'"%s": %d',model2.rxns{i},opt2.x(i));
   if i~=n
       fprintf(file,', ');
   end
end
fprintf(file,'}');
fclose(file);





save('makeEscherInfo.mat');

end

