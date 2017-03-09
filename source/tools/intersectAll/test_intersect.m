close all;
clear;

NnodesPerSet = 1000:100:1000;
NN = numel(NnodesPerSet);

Npow     = 8;
Nsets    = 2^Npow;

% Nsets    = 100;
maxNoIdx = 10000;

timeOld=zeros(Npow,NN);
timeNew=zeros(Npow,NN);

for j=1:NN

    % vecSets=cellfun(@(x) randi(maxNoIdx,1,x),num2cell(randi(maxNoIdx,1,Nsets)),'UniformOutput',false);
    vecSets=cellfun(@(x) randperm(x),num2cell(NnodesPerSet(j)*ones(1,Nsets)),'UniformOutput',false);

    for i=1:Npow
        tOld=tic;
        [indexSet, numAllintersect] = intersectAll(vecSets(1:2^i));
        timeOld(i,j)=toc(tOld);

        tNew=tic;
        [indexSet_my, numAllintersect_my] = intersectAll_my(vecSets(1:2^i));
        timeNew(i,j)=toc(tNew);
    end

    figure;
    plot(2.^(1:Npow),timeOld(:,j),'r',2.^(1:Npow),timeNew(:,j),'g',2.^(1:Npow),timeOld(:,j)./timeNew(:,j),'b');
    title(sprintf('Nodes per Set: %d',NnodesPerSet(j)));
    legend('Old code','New code','Speedup');
end
 

