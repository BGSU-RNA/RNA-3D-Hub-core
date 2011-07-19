% zPredictMutationRate

% File = zAddNTData({'2avy','2aw4','2j01'},0,[],1);
% File = zAddNTData({'2avy','2aw4','2j01'},0,File,1);

% File = zAddNTData({'2j01'},0,[],1);

% File = zAttachAlignment(File,1);

for f = 1:length(File),
 chain = cat(1,File(f).NT.Chain);
 ch = unique(chain);
 for u = 1:length(ch),

  fid = fopen(['MutationPredictionData' File(f).Filename '_' ch(u) '.txt'],'w');
  FastaCol = [];

  for i = 1:length(File(f).NT),
    if isempty(File(f).NT(i).FASTACol),
      FastaCol(i) = 0;
    else
      FastaCol(i) = File(f).NT(i).FASTACol;
    end
  end

  FastaCol = FastaCol .* (chain' == ch(u));

  for c = 1:max(FastaCol),
    i = find(FastaCol == c);
    if isempty(i),
      fprintf(fid,'%5d   .   .   .   .   .   .   .   .\n', c);
    else                                   % there is a NT corresponding to c

      if length(i) > 1,                    % why does this happen?
        i = i(1);
      end

      e = File(f).Edge(i,:);
      b = File(f).BasePhosphate(i,:);      % the base in a BPh interaction
      p = File(f).BasePhosphate(:,i);      % the phosphate in a BPh interaction
      b(i) = 0;
      p(i) = 0;                            % eliminate self BPh interactions
      Data(1) = sum(abs(fix(e)) == 1);    % # cWW basepairs
      Data(2) = sum((abs(fix(e)) > 1) .* (abs(e) < 14)); % # non-cWW basepairs
      Data(3) = sum((abs(fix(e)) > 19) .* (abs(e) < 24)); % # stacking
      Data(4) = sum(b>0);                  % # BPh with i as the base
      Data(5) = sum(p>0);                  % # BPh with i as the phosphate
      Data(6) = 0; % 1 if i is near an amino acid - can't do this!

      fasta = File(f).NT(i).FASTA;
      total = sum(fasta == 'A') + sum(fasta == 'C') + sum(fasta == 'G') + sum(fasta == 'U');



      Data(7) = sum(fasta == File(f).NT(i).Base) / total; % conservation %

      fprintf(fid,'%5d %5d %2d %2d %2d %2d %2d %2d %7.4f\n', c, i, Data(1), Data(2), Data(3), Data(4), Data(5), Data(6), Data(7));


% 1 - column from FASTA file
% 2 - index of NT in 3D structure file
% 3 - #cWW pairs
% 4 - #non-cWW pairs
% 5 - # stacking
% 6 - # BPh with i as the base
% 7 - # BPh with i as the phosphate
% 8 - 1 if near an amino acid; all 0's here!
% 9 - conservation percentage, 0 to 1.

    end
  end

  fclose(fid);
 end
end


