
% if File is not already defined, set File = []; before running this script

if ~exist('File'),
  File = [];
end

fprintf('Nov 1 No 1 model:\n');
[File,Node] = pMotifModel(File,'1s72','1095','1261','1100');
fprintf('\n');



fprintf('Nov 1 No 2 model:\n');
[File,Node] = pMotifModel(File,'1j5e','580','761','585');
fprintf('\n');



fprintf('Nov 1 No 3 model:\n');
%[File,Node] = pMotifModel(File,'1j5e','1303','1334','1308');
[File,Node] = pMotifModel(File,'2aw4','703','728','708');
fprintf('\n');



%fprintf('Nov 1 No 4 model:\n');
%[File,Node] = pMotifModel(File,'1s72','705','723','710');
%fprintf('\n');

