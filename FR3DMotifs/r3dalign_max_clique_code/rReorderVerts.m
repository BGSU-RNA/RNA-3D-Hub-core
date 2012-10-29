%arranges vertices in ascending or descending degree order, depending on
%dir
function [S] = rReorderVerts(Z,dir)
if nargin < 2
    'ERROR - must input argument for ordering direction'
    return;
end

if isequal(dir,'ascend')
   [C S]=sort(full(sum(Z)),2,'ascend');
elseif isequal(dir,'descend')
   [C S]=sort(full(sum(Z)),2,'descend');
end

