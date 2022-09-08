function [outputArg1,outputArg2] = makeHTMLexample(inputArg1,inputArg2)
%MAKEHTMLEXAMPLE この関数の概要をここに記述
%   詳細説明をここに記述

load('e_coli_core.mat');
load('succ_e.mat');
GRPREXchecker(e_coli_core,targetMet,gvalue);
load('GRPREXchecker.mat');
gdel=find(cell2mat(gvalue(:,2))==0);
d=size(gdel,1);
model=e_coli_core;

file=fopen('e_coli_core&succ_e.html','w');
fprintf(file,'<head>\n');
fprintf(file,'<title>MetNetComp</title>\n');
fprintf(file,'</head>\n');
fprintf(file,'<body>\n')
fprintf(file,'Model = e_coli_core<br>\n');
fprintf(file,'Target metabolite = %s<br>\n',targetMet);
fprintf(file,'Gene deletion size = %d<br>\n',d);
fprintf(file,'Gene deletion: \n');
for i=1:d
    fprintf(file,'%s ',model.genes{gdel(i,1)});
end
fprintf(file,'<br><br>\n');
fprintf(file,'When growth rate is maximized,<br>\n')
fprintf(file,'&nbsp Growth Rate = %f<br>\n',GR);
fprintf(file,'&nbsp Minimum Production Rate = %f<br><br>\n',PR);
fprintf(file,'Substrate:<br>\n');
for i=1:num_sub
    s=sprintf('&nbsp %s',substrate{i,1});
    fprintf(file,s);
    fprintf(file,' : %f<br>\n',-substrate{i,2});
end
fprintf(file,'<br>\n');
fprintf(file,'Product:<br>\n');
for i=1:num_sub
    s=sprintf('&nbsp %s',product{i,1});
    fprintf(file,s);
    fprintf(file,' : %f<br>\n',product{i,2});
end
fprintf(file,'<br>');
fprintf(file,'The JSON file for visualization by Escher: ')
fprintf(file,'<a href="flow.json">Download</a><br>\n');
fprintf(file,'Go to <a href="https://escher.github.io/"> Escher site </a><br>\n');
s1=sprintf('Go to <a href="https://escher.github.io/#/app?map=e_coli_core.Core');
s2=sprintf('%20metabolism&tool=Builder&model=e_coli_core"> Escher site </a><br>\n');

fprintf(file,'</body>\n')
fclose(file);

  


save('makeHTMLexample.mat');

end

