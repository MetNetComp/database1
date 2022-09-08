function [outputArg1,outputArg2] = trimGdelEcolicore(inputArg1,inputArg2)
load('e_coli_core.mat');
load('analyzeResult.mat','minFVA');
model=e_coli_core;
m=size(model.mets,1);

for i=1:m
    i
    if (minFVA(i,1)<0.001) && (minFVA(i,2)>=0.001) && (minFVA(i,3)>=0.001)
        s=sprintf('results/%d.mat',i);
        load(s,'gvalue');
        [gvalue, finalGRPR, size1, size2, size3, time]=trimGdel(model,model.mets{i},gvalue,i);
        result(i,:)=horzcat(finalGRPR,size1,size2,size3);
        ttt(i,1)=time;
    end
end



save('trimGdelEcolicore.mat');

end

