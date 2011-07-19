% zMutationRatePredictors



% File = zAddNTData({'2avy','2aw4','2j01'},0,[],1);
% File = zAddNTData({'2avy','2aw4','2j01'},0,File,1);

% File = zAddNTData({'2j01'},0,[],1);

% File = zAttachAlignment(File,1);

% ---------- incorporate Jesse's list of bases near a protein (for some files)

for f = 1:length(File),
  for i = 1:length(File(f).NT),
    File(f).NT(i).BaseProtein = '.';
  end

  ff = find(ismember({'2avy','2aw4'},lower(File(f).Filename)));
  if ~isempty(ff),
    [n,t] = xlsread('Bases_from_Ec_and_Interactions_10_21-1.xls');
    for r = 1:length(n(:,1)),
      if strcmpi(File(f).Filename,t{r,1}),
        switch ff,
          case 1, i = zIndexLookup(File(f),num2str(n(r,1)),'A');
          case 2, i = zIndexLookup(File(f),num2str(n(r,1)),'B');
        end
        if length(i) == 1,
          if n(r,15) == 1,
            File(f).NT(i).BaseProtein = '1';
          else
            File(f).NT(i).BaseProtein = '0';
          end
        end
      end
    end
  end
end

% ----------------------------------------------- Determine center of PDB file

% E coli: (2avy, 2aw4)
moleculecenter = [24.1858  157.7919 -120.8290];

% T.th: (2j01, 2j00, not 1j5e)

% moleculecenter = [-27.9808  118.6233  138.5073];

for f = 1:length(File),

% centers = cat(1,File(f).NT.Center);
% moleculecenter = mean(centers);


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
      fprintf(fid,'%5d .        .     .  .  .  .  .  .  .\n', c);
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

      fasta = File(f).NT(i).FASTA;
      total = sum(fasta == 'A') + sum(fasta == 'C') + sum(fasta == 'G') + sum(fasta == 'U');


%      Data(7) = sum(fasta == File(f).NT(i).Base) / total; % conservation %

      Data(8) = norm(moleculecenter - File(f).NT(i).Center);

      fprintf(fid,'%5d %1s %8s %5d %2d %2d %2d %2d %2d %2c\n', c, File(f).NT(i).Base, File(f).NT(i).Number, i, Data(1), Data(2), Data(3), Data(4), Data(5), File(f).NT(i).BaseProtein);


% 1 - column from FASTA file
% 2 - base letter in 3D structure
% 3 - nucleotide number in 3D structure
% 4 - index of NT in 3D structure file
% 5 - #cWW pairs
% 6 - #non-cWW pairs
% 7 - # stacking partners
% 8 - # BPh with i as the base
% 9 - # BPh with i as the phosphate
%10 - 1 if near an amino acid

    end
  end

  fclose(fid);
 end
end


