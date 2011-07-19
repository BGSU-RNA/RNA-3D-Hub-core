% zFASTAtoRDFAlignment(FASTA) writes out a FASTA-format alignment as an RDF-format alignment, presuming that each column is a NT to NT correspondence, a simplifying assumption

% zReadFASTA reads the indicated fasta file and stores the records in an array.  Each element has three fields:
% Data(n).Header is the header for the nth sequence
% Data(n).Aligned is the sequence, with gaps
% Data(n).Sequence is the sequence with all gaps stripped out

% FASTA = zReadFASTA('Alignments/5S_Bacterial_Stombaugh_et_al_Sup_Mat_S1.fasta')
% Text = zFASTAtoRDFAlignment('5S_Bacterial_Stombaugh_et_al_Sup_Mat_S1.fasta',FASTA);

function [Text,Group] = zFASTAtoRDFAlignment(Filename,FASTA,Verbose)

if nargin < 2,
  FASTA = zReadFASTA(Filename); 
end

if nargin < 3,
  Verbose = 1;
end

for s = 1:length(FASTA);
  if Verbose > 0,
    fprintf('%5d %s\n', s, FASTA(s).Aligned);
  end
  Array(s,:) = FASTA(s).Aligned;              % normal array with gaps
  Array(s,:) = strrep(Array(s,:), '.', '-');  % use - as gap character
end

if any(Array(1,:) == '(') && any(Array(1,:) == ')'),   % dot-bracket notation
  PairColumn = pParseDotBracket(Array(1,:),1);
  startseq = 2;
else
  PairColumn = [];
  startseq = 1;
end

[S,C] = size(Array);

for s = startseq:S,
  ColToPos(s,:) = cumsum(Array(s,:) ~= '-');
end

PrefixCode = 2;

r = 1;
Text{r} = '@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>.';
r = r + 1;
Text{r} = '@prefix seq : <http://database_of_seqences.org/>.'; % Fictitious database of sequences
r = r + 1;
Text{r} = '@prefix alo : <http://purl.obolibrary.org/alo/>.'; % Fictitious location of structural alignment ontology
r = r + 1;
Text{r} = '@prefix rnao: <http://purl.obolibrary.org/obo/>.';
r = r + 1;


AligName = Filename;
AligName = strrep(AligName,'.fasta','.n3');
AligName = strrep(AligName,'.FASTA','.n3');
  AligName = strrep(AligName,' ','_');
  AligName = strrep(AligName,':','_');
  AligName = strrep(AligName,'<','_');
  AligName = strrep(AligName,'>','_');
  AligName = strrep(AligName,'?','_');
  AligName = strrep(AligName,'*','_');
  AligName = strrep(AligName,'&','_');

Text{r} = ['@prefix     : <http://rna.bgsu.edu/Sequence_Alignments/' AligName '/>.'];    % This forms the namespace for the current document
r = r + 1;


% ------------------------------------------ gather prefixes for sequences

for s = startseq:S,                        % go through each sequence
  SN = FASTA(s).Header;
  SN = strrep(SN,' ','_');
  SN = strrep(SN,':','_');
  SN = strrep(SN,'<','_');
  SN = strrep(SN,'>','_');
  SN = strrep(SN,'?','_');
  SN = strrep(SN,'*','_');
  SN = strrep(SN,'&','_');

  SURI{s} = SN;

  SN = sprintf('%6d', s);                 
  SN = strrep(SN,' ','0');                       % pad with leading zeros
  SN = ['Seq' SN];

  SName{s} = SN;

  Text{r} = ['@prefix ' SN ':     <http://rfam.org/RDF/' SURI{s} '>.'];
                                               % prefix for each sequence
  r = r + 1;
end

% ------------------------------------------ set up correspondence groups

g = 0;                                     % group counter

for c = 1:C,                               % go through each column

  if any(Array(:,c) ~= '-'),               % non-gap character in this column
  
    GroupName{c} = [':NTNT' sprintf('%5d',c)];    % name of NT-NT corresp group
    GroupName{c} = strrep(GroupName{c},' ','0');  % pad with leading zeros

    g = g + 1;
    Group(g).Instance = [];                                % initialize this group
    Group(g).Name = GroupName{c};

    Text{r} = [GroupName{c} ' rdfs:type rnao:nucleotide_correspondence_group.'];
    r = r + 1;

    Text{r} = [GroupName{c} ' rdfs:label Nucleotide_nucleotide_correspondence.'];
    r = r + 1;

    for s = 1:S,                               % loop through sequences
      if Array(s,c) ~= '-',                    % if not a gap,
        Text{r} = [SName{s} '/' num2str(ColToPos(s,c)) ' rnao:part_of ' GroupName{c} '.'];
        r = r + 1;
        Group(g).Instance(s) = ColToPos(s,c);
      else
        Group(g).Instance(s) = 0;
      end
    end
  end
end

% ------------------------------------------ 

[aa,bb] = size(PairColumn);

for p = 1:aa,                                  % loop through pairs
  c = PairColumn(p,1);
  d = PairColumn(p,2);
  PairName{p} = [':Pair' sprintf('%5d',c)];    % name of Pair corresp group
  PairName{p} = strrep(PairName{p},' ','0');  % pad with leading zeros

  Text{r} = [PairName{p} ' rdfs:type rnao:cww_basepair.'];
  r = r + 1;

  Text{r} = [PairName{p} ' rdfs:label cWW_basepair.'];
  r = r + 1;
  
  Text{r} = [GroupName{c} ' rnao:first_part_of ' PairName{p} '.'];
  r = r + 1;

  Text{r} = [GroupName{d} ' rnao:second_part_of ' PairName{p} '.'];
  r = r + 1;

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

