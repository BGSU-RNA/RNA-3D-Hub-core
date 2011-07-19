% zAlignToFASTA determines a vector FastaCol which maps the indices in File corresponding to Chain to the corresponding column numbers in the Aligned field of the data in FASTA.  If Entry is specified or non-zero, zAlignToFASTA uses the indicated entry of the FASTA file as the key.  Otherwise, it tries aligning to each sequence successively until a good match is found.

% The presumption is that an alignment contains information from only one chain

% File = zAddNTData('2avy');
% FASTA = zReadFASTA('Alignments\16S_Bacterial_Stombaugh_et_al_Sup_Mat_S2.fasta');
% Entry = 1;

% F = zAlignToFASTA(File,'A',FASTA,28,2);


function [File] = zAlignToFASTA(File,Chain,FASTA,Entry,Verbose)

if nargin < 5,
  Verbose = 0;
end

if nargin < 4,
  Entry = 0;
elseif Entry == 0,
  Verbose = 2;
end

C = cat(2,File.NT.Chain);
a = find(C == Chain);

StructureSequence = cat(2,File.NT(a).Base);

% ------------------------------------ Identify matching alignment, if needed

if Entry == 0,
  if Verbose > 0,
    fprintf('Aligning sequence from structure to each entry in FASTA file\n');
  end
  for e = 1:length(FASTA),
    p = 0.99;
    d = 2;
%   [matches,align1,align2,s1,s2] = dNeedlemanWunsch(StructureSequence,FASTA(e).Sequence,p,d);

    s = FASTA(e).Sequence;
    s = strrep(s,'?','N');

    [matches,align1,align2,s1,s2] = zNeedlemanWunsch(StructureSequence,s);

    if Verbose > 1,
      fprintf('Sequence %4d has %4d matches and characters %8s; from %s\n', e, matches, unique(s), FASTA(e).Header);
    end

    m(e) = matches;
  end

  [y,Entry] = max(m);

  if Verbose > 1,
    [y,i] = sort(-m);
    fprintf('Top choices:\n');
    for j = 1:10,
      fprintf('Sequence %4d has %4d matches and characters %8s; from %s\n', i(j), m(i(j)), unique(FASTA(i(j)).Sequence), FASTA(i(j)).Header);
    end
  end      

  if Verbose > 0,
    fprintf('Using sequence %d, which has %d matches out of %d bases.\n', Entry, m(Entry), length(StructureSequence));
  end

%%%%% Note: this might not be the right file!

  [n,t,r] = xlsread([DropboxRoot filesep 'Alignments\StructureToAlignmentMap.xls']);

  for i = 1:length(r(:,1)),
    if strcmp(upper(File.Filename),upper(r{i,1})) && (Chain == r{i,2}),
      r{i,4} = Entry;
    end
  end


%%%%%% Note:  the file being written here may not be the file you are reading!

  xlswrite([DropboxRoot filesep 'Alignments\StructureToAlignmentMap.xls'],r);

  if Verbose > 0,
    fprintf('Added this choice to structure to alignment map.\n');
  end
end

% ------------------------------------ Align 3D structure sequence to alignment

  % --------- Align 3D structure sequence simply to sequence from alignment

  p = 0.99;
  d = 2;
%  [matches,align1,align2,s1,s2] = dNeedlemanWunsch(StructureSequence,FASTA(Entry).Sequence,p,d);

  [matches,align1,align2,s1,s2] = zNeedlemanWunsch(StructureSequence,FASTA(Entry).Sequence,8,4);

  k = double(FASTA(Entry).Aligned ~= '-');            % 1 if not a gap
  m = cumsum(k);
  
  x = zInvertFunction(align1);                        % map 

  y = [align2 (length(FASTA(Entry).Sequence)+1)]; 
  z = [zInvertFunction(m) length(FASTA(Entry).Aligned)]; 
  % ---------- z maps elements of FASTA sequence to columns of the alignment

%align1

%length(x)
%length(y)
%length(z)

  FastaCol = z(y(x));                     % maps index in 3D structure to
                                          % column in FASTA file

  Alignment = cat(1,FASTA.Aligned);       % put together alignment into matrix

%FastaCol
%min(FastaCol)
%max(FastaCol)

  L = length(Alignment(1,:));                               % max num columns

  for i = a,
    colnum  = FastaCol(min(length(FastaCol),i-a(1)+1));                   % current FASTA col
    nextgap = min(colnum+1,L);                      % next column
    lastgap = min(FastaCol(min(length(FastaCol),i+1-a(1)+1))-1,L);        % 
    File.NT(i).FASTA     = Alignment(:,colnum);     % pull out column
    File.NT(i).FASTACol  = colnum;                  % record column #
    File.NT(i).GapsAfter = Alignment(:,nextgap:lastgap);
                % pull out columns between those aligned to the 3D structure
  end 

  if Verbose > 1,
    fprintf('Structure: %s\n', StructureSequence);
    fprintf('Alignment: %s\n', FASTA(Entry).Aligned(FastaCol));

    clear c d
    for i = 1:length(File.NT),
      if length(File.NT(i).FASTA > 0),
        c(1,i) = File.NT(i).Base;
        d(1,i) = File.NT(i).FASTA(Entry);
      end
    end
    fprintf('\n');
    fprintf('Structure: %s\n', c);
    fprintf('Alignment: %s\n', d);

    % show all columns of Entry from alignment, with nucleotides from 3D struct

    for i = 1:length(FASTA(Entry).Aligned),
      StructAligned(i) = '-';
    end

    for i = 1:length(FastaCol),
      StructAligned(FastaCol(i)) = StructureSequence(i);
    end

    i = find( (StructAligned ~= '-') + (FASTA(Entry).Aligned ~= '-'));

    fprintf('\n');
    fprintf('Structure: %s\n', StructAligned(i));
    fprintf('Alignment: %s\n', FASTA(Entry).Aligned(i));
  end

