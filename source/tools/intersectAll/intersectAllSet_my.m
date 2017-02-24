function [vals] = intersectAllSet_my(sets)
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

    % make cell-array of rows
    sets = cellfun(@(x) x(:)',sets,'UniformOutput',false);

    Ns      = numel(sets);
    setsIdx = num2cell(1:Ns);

    I = cellfun(@(x,y) repmat(y,1,numel(x)),sets,setsIdx,'UniformOutput',false); I = [I{:}];
    J = [sets{:}];
    M = sparse(I,J,true,Ns,max(J));

    vals = sets{1}(all(M(:,sets{1}),1));

%   numSet          = length(sets);
%   if(~ numSet > 1) 
%       error('Set size must be at least 2 sets'); 
%   end
%   vals = [];
%   vals = intersect(sets{1}, sets{2});
%   if(numSet > 2)
%     for i=3:numSet
%       vals = intersect(vals, sets{i});
%     end
%   end
end