% zNTtoURI takes one or more indices and a 3D structure file and outputs a URI

function [URI,Prefix] = zNTtoURI(File,NTList,PrefixCode)

if nargin < 3,
  PrefixCode = 0;
end

if nargin < 4,
  Abbreviation = File.Filename;
end

switch PrefixCode,
case 0,
  Prefix = '';
case 1,
  Prefix = ['@prefix ' Abbreviation ':       <http://pdb.org/RDF/>.'];
case 2,
  Prefix = ['@prefix ' Abbreviation ':       <http://pdb.org/RDF/' File.Filename '/>.'];
end

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

if nargin == 1 || isempty(NTList),
  NTList = 1:length(File.NT);                  % display them all
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(NTList),'char'),
  NTList = {NTList};
end

if strcmp(class(NTList),'cell'),
  Indices = zIndexLookup(File,NTList);
else
  Indices = NTList;
end

% ---------------------------------------------- Create URIs

for i = 1:length(Indices),
  switch PrefixCode,
  case 0,
    URI{i} = ['http://pdb.org/RDF/' File.Filename '/'];
  case 1,
    URI{i} = [Abbreviation ':' File.Filename '/'];
  case 2,
    URI{i} = [Abbreviation ':'];
  end

  j = Indices(i);

  if ~isfield(File.NT,'Model'),
    File.NT(j).Model = 1;                        % hack
  end

  URI{i} = [URI{i} num2str(File.NT(j).Model) '/' File.NT(j).Chain '/' File.NT(j).Number];
end

