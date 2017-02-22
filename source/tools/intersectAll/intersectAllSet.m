function [intersectValues] = intersectAllSet(vecSets)
%
%  INTERSECTALLSET given set of sets, produce intersection over all sets
%               
%
%  INPUT  :
%    vecSets: cell array of vectors
%
%  OUTPUT : 
%    intersectValues : Vector, values of intersection
%
%  Author : Mehmet Suzen <msuzen at gmail com>  
%  License: BSD
% 
%  Example:
%    vecSets          = {[ 1 3 4 5 6 7 9 ],  [ 1 2 3 6 7 8 9 ], [ 2 3 4 5 6 7 8 ] }
%    intersectValues  = intersectAllSet(vecSets);
%    ans =
%    3     6     7
%    
%
  numSet          = length(vecSets);
  if(~ numSet > 1) 
      error('Set size must be at least 2 sets'); 
  end
  intersectValues = [];
  intersectValues = intersect(vecSets{1}, vecSets{2});
  if(numSet > 2)
    for i=3:numSet
      intersectValues = intersect(intersectValues, vecSets{i});
    end
  end
end