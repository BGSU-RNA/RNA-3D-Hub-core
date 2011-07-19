% xSearchtoRDFAlignment(Search) writes out the alignment of 3D motifs in RDF as an alignment between nucleotides in a 3D structure file

function [Text] = xSearchtoRDFAlignment(Search)

PrefixCode = 2;

Cand = Search.Candidates;
[L,N] = size(Cand);
N = N - 1;                                 % number of nucleotides

FN = unique(Cand(:,N+1));                  % unique file numbers

r = 1;
Text{r} = '@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>.';
r = r + 1;
Text{r} = '@prefix seq : <http://database_of_seqences.org/>.'; % Fictitious database of sequences
r = r + 1;
Text{r} = '@prefix alo : <http://purl.obolibrary.org/alo/>.'; % Fictitious location of structural alignment ontology
r = r + 1;
Text{r} = '@prefix rnao: <http://purl.obolibrary.org/obo/>.';
r = r + 1;

if isfield(Search,'Name'),
  AligName = Search.Name;
else
  AligName = ['example_alignment_' datestr(now,31) '.n3'];
  AligName = strrep(AligName,' ','_');
  AligName = strrep(AligName,':','_');
  AligName = strrep(AligName,'<','_');
  AligName = strrep(AligName,'>','_');
  AligName = strrep(AligName,'?','_');
  AligName = strrep(AligName,'*','_');
  AligName = strrep(AligName,'&','_');
end

Text{r} = ['@prefix     : <http://rna.bgsu.edu/3D_Motif_Alignments/' AligName '/>.'];    % This forms the namespace for the current (fictitious) document
r = r + 1;


% ------------------------------------------ gather prefixes for files

for f = 1:length(FN),
  [URI,Prefix] = zNTtoURI(Search.File(FN(f)),1,PrefixCode);
  Text{r} = Prefix;                        % prefix for each distinct file
  r = r + 1;
end

% ------------------------------------------ set up correspondence groups

for n = 1:N,
  GroupName = [':NTNT' sprintf('%4d',n)];       % name of NT-NT corresp group
  GroupName = strrep(GroupName,' ','0');           % pad with leading zeros

  Text{r} = [GroupName ' rdfs:type rnao:nucleotide_correspondence_group.'];
  r = r + 1;

  Text{r} = [GroupName ' rdfs:label Nucleotide_nucleotide_correspondence.'];
  r = r + 1;

  for c = 1:L,                               % loop through candidates
    f = Cand(c,N+1);                         % file number
    T = zNTtoURI(Search.File(f),Cand(c,n),PrefixCode);  % URI of this nucleotide
    Text{r} = [T{1} ' rnao:part_of ' GroupName '.'];
    r = r + 1;
  end
end

% ------------------------------------------ 

for r = 1:length(Text),
  fprintf('%s\n', Text{r});
end

% ------------------------------------------ 

fid = fopen(AligName,'w');
for r = 1:length(Text),
  fprintf(fid,'%s\n', Text{r});
end
fclose(fid);

