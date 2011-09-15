% pReadFASTA(Filename) reads a Filename.fasta and organizes the sequences
% and organism names into a vector Sequence

% The user may specify a Key Organism and a Key Sequence, which makes it
% easier to know which columns of the FASTA file are desired

function [Sequence] = pReadFASTA(Filename,A,B)

if nargin == 1,                         % no A, B specified
  A = 1;
  B = Inf;                              % fix this later
end

if strcmp(class(A),'char'),
  KeyOrg = A;
  KeySeq = B;
else
  First = A;                            % first column in FASTA file
  Last  = B;                            % last column
  KeyOrg = 'Z';                         % just a dummy
end

fid = fopen(Filename,'r');

if fid > 0,

i = 0;
L = 1;

while L > -1
  L = fgetl(fid);
  if L > -1
    s = strfind(L, '>');
    if ~isempty(s)
      i = i + 1;
      Sequence(i).Organism = L((s+1):(length(L)));
      Sequence(i).Fasta    = '';
    else
      Sequence(i).Fasta = [Sequence(i).Fasta L(1:(length(L)))];
    end
  end
end

fclose(fid);

KeyIndex = 0;

for k=1:length(Sequence),                           % loop through sequences
  Sequence(k).Fasta = strrep(Sequence(k).Fasta,'.','-');

  s = Sequence(k).Fasta;                            % fasta sequence w/ gaps
  [m,i] = sort((s=='-')+(1:length(s))/(2*length(s))); % i lists non-gap columns
  L = sum(s ~= '-');                                % num of non-gap charact.
  Sequence(k).FastaIndex = i(1:L);                  
          % map from ungapped sequence positions to column number in Fasta file
  Sequence(k).IndextoColumn = i(1:L);                  
          % map from ungapped sequence positions to column number in Fasta file
  Sequence(k).ColumntoIndex = max(1,cumsum(s ~= '-'));
          % map from column number in Fasta file to ungapped sequence position
  Sequence(k).FastaNum = num2str(k);                % entry number in file
  if ~isempty(strfind(Sequence(k).Organism,KeyOrg)),
    KeyIndex = k;
  end
end


if (KeyIndex > 0) & strcmp(class(A),'char'),
  a = strfind(Sequence(KeyIndex).Fasta(Sequence(KeyIndex).FastaIndex),KeySeq);

  if ~isempty(a),

    b = a + length(KeySeq) - 1;

    First = Sequence(KeyIndex).FastaIndex(a);          % starting column
    Last  = Sequence(KeyIndex).FastaIndex(b);          % ending column

  else
    fprintf('Could not find key sequence in %s.fasta\n', Filename);
  end

elseif strcmp(class(A),'char'),

  fprintf('Could not find key organism in %s.fasta\n', Filename);

end

if Last == Inf,
  Last = length(Sequence(1).Fasta);
  a    = 1;
else
  a = 1;
end

if ((KeyIndex > 0) & ~isempty(a)) | ~strcmp(class(A),'char'),
  for k = 1:length(Sequence),
    x = Sequence(k).Fasta;
    x = x(First:Last);                          % extract relevant columns
    x = strrep(x,'-','');                       % discard insertion markers
    Sequence(k).X = x;                          % save shortened sequence
  end

  fprintf('Read %s\n', Filename)
end

else

  fprintf('pReadFASTA could not open file %s\n', Filename);

end