function [indexSet, numAllintersect] = intersectAll(vecSets)
%
%  INTERSECTALL given set of sets, produce all possible intersections 
%               for each set.
%
%  INPUT  :
%    vecSets: cell array of vectors
%
%  OUTPUT : 
%    indexSet: cell array of vectors 
%              (each indices of a set that have
%              intersections to all other sets.
%    numAllintersect: 
%              vector number of intersections to all other sets for
%              each set
%  Author : Mehmet Suzen <msuzen at gmail com>  
%  License: BSD
% 
%  Example:
%    vecSets = {[11, 12, 22, 14], [11, 15, 17, 22, 33], [ 17, 25, 14] }
%    [indexSet, numAllintersect] = intersectAll(vecSets);
%    indexSet{1}
%    ans =
%       1     3     4
%    indexSet{2}
%    ans =
%     1     3     4
%      indexSet{3}
%    ans =
%     1     3
%    numAllintersect    
%     3
%     3
%     2
%
  numSet    = length(vecSets);
  allPairs  = nchoosek(1:numSet, 2);
  indexSet  = cell(numSet, 1);
  for i=1:length(allPairs(:,1))
    pairI   = allPairs(i, :);
    iFirst  = vecSets{pairI(1)};
    iSecond = vecSets{pairI(2)};
    [~, iOne, iTwo] = intersect(iFirst, iSecond);
    indexSet{pairI(1),pairI(2)} = [iOne];
    indexSet{pairI(2),pairI(1)} = [iTwo];
%     indexSet{pairI(1),i} = [indexSet{pairI(1)}', iOne];
%     indexSet{pairI(2),i} = [indexSet{pairI(2)}', iTwo];
  end
   % indexSet        = cellfun(@(vv) unique(vv), indexSet, 'UniformOutput', false);
    numAllintersect = cellfun(@(vv) length(vv), indexSet);
end
