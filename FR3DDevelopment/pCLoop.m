
% if File is not already defined, set File = []; before running this script

if ~exist('File'),
  File = [];
end

fprintf('C-Loop model:\n');
[File,Node] = pMotifModel(File,'2aw4','2680','2727','2685');
fprintf('\n');

