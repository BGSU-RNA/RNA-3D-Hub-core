% pMSASequencesJAR3D reads one fasta file corresponding to a motif, eliminates sequences found in 3D structures, parses a fixed number of sequences from it against all models, ranks them, and displays diagnostics.  It uses 2, 3, 4, ... randomly chosen sequences and displays diagnostics of the rank of the correct model.

% Example:
% Scores = JAR3DMatlab.MotifParseSingle(pwd,'IL_225_test2.fasta','IL_225_05_cWW-cWW-cSH');
% S = JAR3DMatlab.MotifTestGeneral(pwd,'IL','IL_test.txt','IL_test.txt');

SeqLimit = 1000;                               % maximum number of seqs to parse

disp('Make sure the Matlab current folder has a MotifLibrary in it');

if ~exist('loopType'),
  disp('Please specify a loop type, for example loopType = ''IL'';')
  break
end

if ~exist('SequenceSource'),
  SequenceSource = 2;                          % parse separate sequences?
end

switch SequenceSource,
case 0,
  sequenceNameFile = [loopType '_Sequences.txt'];
case 1,
  sequenceNameFile = [loopType '_SeparatedSequences.txt'];
  disp('Make sure that you have run pSeparateSearches');
case 2,
  sequenceNameFile = [loopType '_SequencesMSA.txt'];
end

% ---------------------------------------- Read model names

clear Names

n = 1;
fid = fopen(['Models' filesep loopType '_Models.txt'],'r');
if fid > 0,
  L = 1;
  while L > -1,
    L = fgetl(fid);                      % read a line
    if L > -1,
      Names{n} = L;
      switch loopType
      case 'IL'
        Names{n+1} = [L ' reversed'];    % sequences have been parsed twice
        n = n + 2;
      case 'HL'
        n = n + 1;
      end
    end
  end
  fclose(fid);
else
  fprintf('Could not open file of model names\n');
  break
end

switch loopType
case 'IL'
  t = length(Names);
  ModelNames = Names(1:2:(t-1));
case 'HL'
  ModelNames = Names;
end

% ---------------------------------------- Read sequence names

n = 1;

fid = fopen(['Sequences' filesep sequenceNameFile],'r');

if fid > 0,
  L = 1;
  while L > -1,
    L = fgetl(fid);                      % read a line
    if L > -1,
      SequenceNames{n} = L;
      n = n + 1;
    end
  end
  fclose(fid);
else
  fprintf('Could not open file of sequence names, using model names instead\n');
end

% ----------------------------------------- Cycle through sequence files

for f = 1:length(SequenceNames),
  FN = SequenceNames{f};                  % current fasta filename

  FASTA = zReadFASTA(['Sequences' filesep FN]);

  % --------------------------------------- write out smaller FASTA file

  canon = {'AU','UA','GC','CG','GU','UG'};   % canonical pairs

  if strcmp(loopType,'IL'),
    fid = fopen(['Sequences' filesep 'MSAtemp.fasta'],'w');
    fidr = fopen(['Sequences' filesep 'MSAtemp_reversed.fasta'],'w');
    n = 1;
    NumSeq = 0;
    while n <= min(length(FASTA),SeqLimit),
      tem = upper(FASTA(n).Aligned);         % sequence with gaps
      i = strfind(tem,'*');
%      if ismember(tem([1 end]),canon) && ismember(tem([i-1 i+1]),canon),
      if ismember(tem([1 end]),canon) && ismember(tem([i-1 i+1]),canon),
        fprintf(fid,'> %s\n', FASTA(n).Header);
        fprintf(fid,'%s\n', tem);
        fprintf(fidr,'> %s\n', FASTA(n).Header);
        t = [tem((i+1):end) '*' tem(1:(i-1))];
        fprintf(fidr,'%s\n', t);
        NumSeq = NumSeq + 1;
      end
      n = n + 1;
    end
    fclose(fid);
    fclose(fidr);
  else
    fid = fopen(['Sequences' filesep 'MSAtemp.fasta'],'w');
    for n = 1:min(length(FASTA),SeqLimit),
      fprintf(fid,'> %s\n', FASTA(n).Header);
      fprintf(fid,'%s\n', FASTA(n).Sequence);
    end
    fclose(fid);
    NumSeq = n;
  end

  SFN  = ['MSAtemp.fasta'];
  SFNR = ['MSAtemp_reversed.fasta'];

  JAR3D_path

  S = [];

  if NumSeq > 0,
    for m = 1:length(ModelNames),
      MFN = ModelNames{m};
      Scores = JAR3DMatlab.MotifParseSingle(pwd,SFN,MFN);
      S = [S Scores];
      if strcmp(loopType,'IL'),
        Scores = JAR3DMatlab.MotifParseSingle(pwd,SFNR,MFN);
        S = [S Scores];
      end
    end

    [y,i] = max(S,[],2);                      % value and location of the max

    [b,t] = zUniqueRows(i);

    MotifNum = FN(4:6);

    for m = 1:length(ModelNames),
      tem = ModelNames{m};
      MN(m) = str2num(tem(4:6));              % model numbers as numbers
      GroupNum{m} = tem(4:6);                 % 3-digit group numbers as strings
      GNtoIndex(MN(m)) = m;
    end

    N = find(ismember(GroupNum,MotifNum));    % number of own model  

    if strcmp(loopType,'IL'),
      NN = 2*N-1;
    else
      NN = N;
    end

    Better = (S > S(:,NN)*ones(1,length(S(1,:))));

    fprintf('Sequence file is %s\n', FN);
    fprintf('Correct model is     %s\n', ModelNames{N});
    fprintf('%d sequences have canonical flanking cWW pairs\n', NumSeq);
    fprintf('Number of times own model has top rank:  %d out of %d\n', sum(i==NN), NumSeq);

    fprintf('Top models are:\n');
    for k = 1:length(t),
      fprintf('%3d times it is:  %s\n', t(k), Names{b(k)});
    end

    fprintf('\n');

    figure(1)
    clf
    spy(Better)
    title('Models which score better than the sequences own model');

    figure(2)
    clf
    maxrank = 20;
    n = hist(min(sum(Better,2),maxrank),1:maxrank);
    bar(1:20,n/sum(n));
    axis([0.5 maxrank+0.5 0 max(n/sum(n))*1.1]);

    set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 12 14 16 18 20])
    set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10','12','14','16','18','>=20'})


    MN = ModelNames{N};
    MN = strrep(MN,'_',' ');
    title([MN ' rank of correct model']);
    MN = ModelNames{N};
    MN = strrep(MN,'.txt','_MSA_hist.png');
    saveas(gcf,MN,'png');

    drawnow

%    pause

  else

    fprintf('No sequences flanked by canonical pairs found for %s\n', FN);

  end

end
