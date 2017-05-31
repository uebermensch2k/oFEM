function [idx, numInt] = intersectAll_my(sets)
%intersectAll computes all intersections for any given set with any other
%
%  [idx, numInt] = intersectAll(sets) given a cell array of vectors,
%  returns all intersections of any set in sets with any other set in sets
%  as idx and the number of elements that actually intersect in numInt.
%  More precisely, provided Ns=numel(sets), idx is a NsxNs cell-array and
%  numInt a NsxNs uin32-array with the following structure:
%
%   - idx{i,j} returns an index array into sets{i} of all elements of
%     intersect(sets{i},sets{j}).
%
%   - numInt{i,j} returns numel(intersect(sets{i},sets{j}))
%
%  It is assumed that the elements in any set in sets are unique but not
%  necessarily sorted, i.e., all(sets{i}==unique(sets{i},'stable'))==1 is
%  true.
%
%  Author : Michael Dudzinski
%  License: BSD
% 
%  Example:
%    sets = {[11, 12, 22, 14], [11, 15, 17, 22, 33], [17, 25, 14]}
%    [idx, numInt] = intersectAll(sets);
%    >> idx{1,:}
%    ans =
%       1   1   1   1
%    ans =
%       1   0   1   0
%    ans =
%       0   0   0   1
%
%     >> idx{2,:}
%     ans =
%        1   0  0   1   0
%     ans =
%        1   1   1   1   1
%     ans =
%        0   0   1   0   0
%
%     >> idx{3,:}
%     ans =
%        0   0   1
%     ans =
%        1   0   0
%     ans =
%        1   1   1
%
%     >> numInt
%     numInt =
%        4   2   1
%        2   5   1
%        1   1   3
%
%     >> sets{1}(idx{1,2})
%     ans =
%         11    22
%
%     >> sets{2}(idx{2,1})
%     ans =
%         11    22
%


%     Ns = numel(vecSets);
% 
%     indexSet = cell(Ns,1);
%     for i=1:Ns
%         idx = [1:i-1,i+1:Ns];
%         indexSet{i}=ismember(vecSets{i},[vecSets{idx}]);
%     end
%     numAllintersect = cellfun(@(vv) length(vv), indexSet);

    % make cell-array of rows
    sets = cellfun(@(x) x(:)',sets,'UniformOutput',false);

    Ns      = numel(sets);
    setsIdx = num2cell(1:Ns);

    I = cellfun(@(x,y) repmat(y,1,numel(x)),sets,setsIdx,'UniformOutput',false); I = [I{:}];
    J = [sets{:}];
    M = sparse(I,J,true,Ns,max(J));
%     M = sparse(I,J,1,Ns,max(J)); % if sets{i} have repeated elements use this

    idx = cellfun(@(x) mat2cell(full(M(:,x)),ones(Ns,1),numel(x)),sets,'UniformOutput',false);
    idx = [idx{:}]';

    numInt = cellfun(@(vv) sum(uint32(vv),'native'), idx);
end
