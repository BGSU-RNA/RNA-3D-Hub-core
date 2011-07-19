
% zReadCraigAlignment reads 3D alignment Excel spreadsheets from Jesse
% Stombaugh.  When that spreadsheet is reformatted, things will need to be
% changed below.

% It should be revised to allow simultaneous alignment of multiple molecules

function [Alignment,File] = zReadCraigAlignment(f1,f2,XLSName,File)

if nargin < 3,
  XLSName = '5S_Tt_Hm_alignment.xls';
end

[n,t,r] = xlsread(XLSName);
fprintf('Read %s\n', XLSName);

Filenames{1} = t{1,1};
Filenames{2} = t{1,5};

for i = 1:length(n(:,1)),
  A{1,i} = {num2str(n(i,1)), num2str(n(i,2))};        % basepairs from file 1
  A{2,i} = {num2str(n(i,5)), num2str(n(i,6))};        % basepairs from file 2
end

Chain{1} = '(B)';
Chain{2} = '(9)';

Filenames = Filenames([f1 f2]);
A         = A([f1 f2],:);
Chain     = Chain([f1 f2]);

if nargin < 4,                               % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(Filenames,2);   % load PDB data
else
  [File,SIndex] = zAddNTData(Filenames,2,File); % add PDB data if needed
end

File = File(SIndex([1 2]));               % pull out the two we need

c = 0;                                      % counter for elements of alignment

clear Alignment

Alignment.Filename{1} = File(1).Filename;
Alignment.Filename{2} = File(2).Filename;

for row = 1:length(A(1,:)),
  if ~strcmp(A{1,row}{1},'NaN') && ~strcmp(A{2,row}{1},'NaN'),
    c = c + 1;
    i1 = zIndexLookup(File(1),[A{1,row}{1} Chain{1}]);
    i2 = zIndexLookup(File(2),[A{2,row}{1} Chain{2}]);
    if ~isempty(i1) && ~isempty(i2),
      Alignment.Correspondence(c).File1(1) = i1;
      Alignment.Correspondence(c).File2(1) = i2;
    end
  end

  if ~strcmp(A{1,row}{2},'NaN') && ~strcmp(A{2,row}{2},'NaN'),
    i1 = zIndexLookup(File(1),[A{1,row}{2} Chain{1}]);
    i2 = zIndexLookup(File(2),[A{2,row}{2} Chain{2}]);
    if ~isempty(i1) && ~isempty(i2),
      Alignment.Correspondence(c).File1(2) = i1;
      Alignment.Correspondence(c).File2(2) = i2;
    end
  end

end

