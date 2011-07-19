% zCompareAlignmentToStructures takes as input a FASTA-format alignment of two (or more?) RNA sequences, together with two (or more?) RNA 3D structures.  It aligns the sequences to the sequences in the 3D structures, then displays bar diagrams illustrating how the aligned nucleotides superimpose in 3D.

function [void] = zCompareAlignmentToStructures(File,FASTA,AlignmentName)

if nargin < 3,
  AlignmentName = 'Unknown';
end

if length(FASTA(:,1)) < 2,
  % ---------------- Treat FASTA as a file name and read

  Data = zReadFASTA(FASTA);

  if nargin < 3,
    i = strfind(FASTA,'.');
    if ~isempty(i),
      AlignmentName = FASTA(1:(i-1));
    else
      AlignmentName = FASTA;
    end
  end

  clear FASTA

  for d = 1:length(Data),
    FASTA(d,:) = Data(d).Aligned;

Data(d).Header

  end
end

Verbose = 1;

NS = length(File);                            % number of sequences/structures

% --------------------------------- Align each line of FASTA file to 3D struct

for f = 1:length(File),
  FASTA(f,:) = upper(FASTA(f,:));
  FASTA(f,:) = strrep(FASTA(f,:),'T','U');

  sequence = FASTA(f,:);
  sequence = strrep(sequence,'-','');         % remove gaps

  structure = cat(2,File(f).NT.Base);

  [matches,align1,align2,s1,s2] = zNeedlemanWunsch(sequence,structure);
  SeqStruct{f}.align1 = align1;
  SeqStruct{f}.align2 = align2;
  SeqStruct{f}.s1 = s1;
  SeqStruct{f}.s2 = s2;

  if Verbose > 0,
    fprintf('Alignment between sequence %d and sequence from %s structure:\n', f, File(f).Filename);
    fprintf('%s\n', s1);
    fprintf('%s\n\n', s2);
  end
end

% ------------------------------- Use FASTA alignment to infer alignment
% ------------------------------- of nucleotides in the 3D structures, pairwise

NTList1 = [];
NTList2 = [];

for i = 1:(NS-1),
  for j = (i+1):NS,

    if Verbose > 0,
      fprintf('Alignment between sequences %d and %d from FASTA:\n', i, j);
      fprintf('%s\n', FASTA(i,:));
      fprintf('%s\n\n', FASTA(j,:));
    end

    ColToSeqiPos = cumsum(FASTA(i,:) ~= '-');    % entry n is position in seq i
    ColToSeqjPos = cumsum(FASTA(j,:) ~= '-');
    for n = 1:length(SeqStruct{i}.align1),% go through seq-struct corresp for i
      StructureiIndex = SeqStruct{i}.align2(n);
      SequenceiPosition = SeqStruct{i}.align1(n);
      k = find(ColToSeqiPos == SequenceiPosition);   % appears in FASTA alignment

      if length(k) > 0,
        FASTAColumn = k(1);

% [n StructureiIndex SequenceiPosition k(1) FASTAColumn]

        if FASTA(j,FASTAColumn) ~= '-',           % aligned to something
          SequencejPosition = ColToSeqjPos(FASTAColumn);
          m = find(SeqStruct{j}.align1 == SequencejPosition);
          if length(m) > 0,
            StructurejIndex = SeqStruct{j}.align2(m(1));

%[n StructureiIndex SequenceiPosition SequencejPosition StructurejIndex 0]

%pause
            NTList1 = [NTList1 StructureiIndex];
            NTList2 = [NTList2 StructurejIndex];
          end
        end
      end
    end

    s1 = cat(2,File(i).NT(NTList1).Base);
    s2 = cat(2,File(j).NT(NTList2).Base);


    if Verbose > 0,
      fprintf('Aligned bases between %s and %s:\n', File(i).Filename, File(j).Filename);
      fprintf('%s\n', s1);
      fprintf('%s\n\n', s2);
    end

    clf
    rBarDiagram(File(i),1:length(File(i).NT),File(j),1:length(File(j).NT),NTList1,NTList2,AlignmentName);

    drawnow

%    pause

    saveas(gcf,[AlignmentName '.pdf']);

  end
end
