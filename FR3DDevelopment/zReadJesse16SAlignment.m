
% zReadJesseAlignment reads 3D alignment Excel spreadsheets from Jesse
% Stombaugh.  When that spreadsheet is reformatted, things will need to be
% changed below.

% It should be revised to allow simultaneous alignment of multiple molecules

function [Alignment,File] = zReadJesseAlignment(f1,f2,XLSName,File)

if nargin < 3,
  XLSName = '3DAlignment\16S_Ec_Tt_Struct_alignment_2_8_07.xls';
end

% do this once:
if ~exist('n') || ~exist('t') || ~exist('r'),
  [n,t,r] = xlsread([pwd filesep XLSName]);
  fprintf('Read %s\n', XLSName);
end

%t{2,11:20}

% These locations are specific to Jesse's 16S alignment spreadsheet

Filenames{1} = t{2,2};
Filenames{2} = t{2,14};
Filenames{3} = t{2,26};

for i = 1:length(n(:,1)),
  A{1,i} = {num2str(n(i, 1)), num2str(n(i, 3))};        % basepairs from file 1
  A{2,i} = {num2str(n(i,13)), num2str(n(i,15))};        % basepairs from file 2
  A{3,i} = {num2str(n(i,25)), num2str(n(i,27))};        % basepairs from file 3
end

Chain{1} = '(A)';
Chain{2} = '(A)';
Chain{3} = '(0)';

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

