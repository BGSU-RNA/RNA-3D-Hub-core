% zFindExemplars finds the best representative for each category of pairs.
% For basepairing, 

NRList = 'Nonredundant_4A_2010-05-19_list';
NRList = 'Nonredundant_4A_2011-01-07_list';

Verbose = 1;               % Verbose = 1 tells it to show distance graphs
LMax    = 500;             % maximum number of pairs to consider in each class

% Pair codes:  1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU

pcodes = [6 7 13 14 15];
pcodes = [6 14 15 16];
pcodes = [6 7 9 11 13 14 15 16];    % pair codes to work on
pcodes = [1 5];
pcodes = [13 14 15 16 1 5 6 7 9 11];    % pair codes to work on

load('PairExemplars','Exemplar');       % load previous exemplars

% load('PairExemplars_Old','Exemplar');   % load last established 

OldExemplar = Exemplar;                 % for when no instances are found

clear Exemplar                          % start fresh

CL = zClassLimits;                      % load basepair classification limits

Pairs{1}  = 'AA';
Pairs{5}  = 'AC';
Pairs{6}  = 'CC';
Pairs{7}  = 'GC';
Pairs{9}  = 'AG';
Pairs{11} = 'GG';
Pairs{13} = 'AU';
Pairs{14} = 'CU';
Pairs{15} = 'GU';
Pairs{16} = 'UU';

% ------------------------------------------ Load non-redundant dataset

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(NRList,0,[],1);   % load PDB data
  File = File(SIndex);
else
  [File,SIndex] = zAddNTData(NRList,0,File,1); % add PDB data if needed
  File = File(SIndex);
end                       

for f = 1:length(File),
  if ~isempty(File(f).Info.Resolution),
    Res(f) = File(f).Info.Resolution;
  else
    Res(f) = 10;
  end

  S(f) = length(File(f).NT);
end

[y,i] = sort(Res);                             % sort files by resolution
File = File(i);
S = S(i);

% ------------------------------------------ Load modeled basepairs

if ~exist('ModelFile'),                        % if no molecule data is loaded,
  [ModelFile,SIndex] = zAddNTData('Model_list',0,[],1);   % load PDB data
else
  [ModelFile,SIndex] = zAddNTData('Model_list',0,ModelFile,1);
end                       

ModelFile = ModelFile(SIndex);

% ------------------------------------------ Load curated exemplars

if ~exist('CuratedFile'),                      % if no molecule data is loaded,
  [CuratedFile,SIndex] = zAddNTData('Curated_list',1,[],1);   % load PDB data
else
  [CuratedFile,SIndex] = zAddNTData('Curated_list',1,CuratedFile,1);
end                       

CuratedFile = CuratedFile(SIndex);

% ------------------------------------------ Make a directory for exemplars

if ~(exist('Exemplars') == 7),        % if directory doesn't yet exist
  mkdir('Exemplars');
end

% loop through paircodes and computer classifications ----------------------

for j = 1:length(pcodes),            % run through all pair codes specified
 pc = pcodes(j);

 CLE = CL(:,1,pc);                   % class limits for this paircode
 CLE = CLE(find(CLE));               % leave out empty entries
 CLE = CLE(find(abs(CLE) < 20));     % leave out stacking

 for row = 1:length(CLE),            % run through classes for this paircode

  % specify criteria for selection of pairs ----------------------------------

  Param = [];
  Param(1,1) = CLE(row);
  Param(1,2) = pc;

  if abs(CLE(row)) < 20,
    Decimal = 1;                          % treat subcategories separately
  else
    Decimal = 0;                          % group subcategories together
  end

  % select pairs using selection criteria -----------------------------------

  List = zFindPairs(CuratedFile,Param,Decimal); % check curated basepairs

  if length(List) > 0,                    % this is a curated basepair

    fprintf('Using a curated pair for %2s %4s class %5.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));

%xDisplayCandidates(CuratedFile,List);

    f  = List(1,3);                       % file number
    Exemplar(row,pc).Filename   = CuratedFile(f).Filename;
    Exemplar(row,pc).Class      = CLE(row);
    Exemplar(row,pc).NT1        = CuratedFile(f).NT(List(1,1));
    Exemplar(row,pc).NT2        = CuratedFile(f).NT(List(1,2));
    Exemplar(row,pc).Resolution = 0;

    CList = zFindPairs(File,Param,Decimal);   % search the non-redundant list

    Exemplar(row,pc).Count      = length(CList(:,1));

    List = zFindPairs(File,Param,Decimal);   % search the non-redundant list

    if length(List) > 0,               % instances of this category are found
      fprintf('FYI, found %5d instances of %2s %4s class %5.1f\n', length(List), Pairs{pc}, zEdgeText(CLE(row),Decimal,pc), CLE(row));

      [PD,ID,i] = zOrderPairs(File,List,LMax,Verbose);

    if 0 > 1,
      if (length(List(:,1)) > 2) && (Verbose > 0),
        zFindExemplarsPlot
      end
      figure(2)
      close
      figure(2)
      VP.Sugar = 1;
      VP.AtOrigin = 1;
      VP.LabelBases = 10;
      zPlotOneNT(Exemplar(row,pc).NT1,VP);
      hold on
      zPlotOneNT(Exemplar(row,pc).NT2,VP);
      axis equal
      rotate3d on

      figure(1)
      close
      figure(1)
      xDisplayCandidates(File,List(i,:));
      rotate3d on
    end

      NT1 = Exemplar(row,pc).NT1;
      NT2 = Exemplar(row,pc).NT2;
      f   = List(i(1),3);
      NT3 = File(f).NT(List(i(1),1));
      NT4 = File(f).NT(List(i(1),2));
      d   = zIsoDiscrepancy(NT1,NT2,NT3,NT4);
      if Verbose > 0,
        fprintf('    IsoDiscrepancy between curated and best from search is %7.4f\n',d);
      end
    end

  else                                      % not a curated basepair

   List = zFindPairs(File,Param,Decimal);   % search the non-redundant list

   if length(List) > 0,               % instances of this category are found

    fprintf('Found %5d instances of %2s %4s class %5.1f ', length(List), Pairs{pc}, zEdgeText(CLE(row),Decimal,pc), CLE(row));

    [PD,ID,i,k] = zOrderPairs(File,List,LMax,Verbose); % order by centrality

    if (length(List(:,1)) > 2) && (Verbose > 0),
      zFindExemplarsPlot                          % display mutual discreps
    end

    f = List(i(1),3);
    Exemplar(row,pc).Filename   = File(f).Filename;
    Exemplar(row,pc).Class      = CLE(row);
    Exemplar(row,pc).NT1        = File(f).NT(List(i(1),1));
    Exemplar(row,pc).NT2        = File(f).NT(List(i(1),2));
    Exemplar(row,pc).Count      = length(List(:,1));
    Exemplar(row,pc).Resolution = File(f).Info.Resolution;

    if 0 > 1,
      figure(1)
      close
      figure(1)
      clf
      xDisplayCandidates(File,List(i,:));
      rotate3d on
    end

   else                                           % no instance found

      List = zFindPairs(ModelFile,Param,Decimal); % check modeled basepairs

      if length(List) > 0,                        % model was found

        fprintf('Using a model for       %2s %4s class %5.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));

        f = List(1,3);
        Exemplar(row,pc).Filename   = ModelFile(f).Filename;
        Exemplar(row,pc).Class      = CLE(row);
        Exemplar(row,pc).NT1        = ModelFile(f).NT(List(1,1));
        Exemplar(row,pc).NT2        = ModelFile(f).NT(List(1,2));
        Exemplar(row,pc).Count      = 0;
        Exemplar(row,pc).Resolution = 0;
      else
        fprintf('No instances and no model for %2s %4s %6.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));

%        pc = 1:16;
        Code1 = mod(pc-1,4)+1;
        Code2 = floor((pc-1)/4)+1;
%        [pc' Code1' Code2']


        [N1,N2,E] = zGetExemplar(CLE(row),Code1,Code2,OldExemplar);

        if ~isempty(E.Filename),
          Exemplar(row,pc).Filename   = E.Filename;
          Exemplar(row,pc).Class      = E.Class;
          Exemplar(row,pc).NT1        = E.NT1;
          Exemplar(row,pc).NT2        = E.NT2;
          Exemplar(row,pc).Count      = E.Count;
          Exemplar(row,pc).Resolution = NaN;

          fprintf('Using previous exemplar for this class\n');
E
        else
          fprintf('***** No previous exemplar found, time to make a new one.\n')
        end

      end
    end
  end

  % add information to speed up the discrepancy calculation later

  [s,t] = size(Exemplar);
  if (row <= s) && (pc <= t),
   if ~isempty(Exemplar(row,pc).NT1),
    E = Exemplar(row,pc);
    Exemplar(row,pc).R            = E.NT2.Rot' * E.NT1.Rot;
    Exemplar(row,pc).T1           = (E.NT2.Center - E.NT1.Center) * E.NT1.Rot;
    Exemplar(row,pc).T2           = (E.NT1.Center - E.NT2.Center) * E.NT2.Rot;
    Exemplar(row,pc).AngleWeight  = [1 1];
    Exemplar(row,pc).LDiscCutoff  = Inf;
   end
  end 

  if LMax >= 500,
    save(['FR3DSource' filesep 'PairExemplars'],'Exemplar'); % Matlab version 7 only
    save PairExemplars_Version_6.mat Exemplar -V6 % for compatibility with older versions
  end

 end
end

% ----------------------------------------- Fix model fields

[s,t] = size(Exemplar);

clear NewExemplar

for i = 1:s,
  for j = 1:t,
    if ~isempty(Exemplar(i,j).Filename),
      E = Exemplar(i,j);
      if isfield(E.NT1,'ModelNum'),
        E.NT1.Model = E.NT1.ModelNum;
        E.NT2.Model = E.NT2.ModelNum;
        E.NT1 = rmfield(E.NT1,'ModelNum');
        E.NT2 = rmfield(E.NT2,'ModelNum');
      end
      if ~isfield(E.NT1,'Model'),
        E.NT1.Model = 1;
        E.NT2.Model = 1;
      end
      if ~isfield(E.NT1,'Unit'),
        E.NT1.Unit = E.NT1.Base;
        E.NT2.Unit = E.NT2.Base;
      end
      E.NT1 = orderfields(E.NT1);
      E.NT2 = orderfields(E.NT2);
      E = orderfields(E);
      NewExemplar(i,j) = E;
    end
  end
end

Exemplar = NewExemplar;

zExemplarIDICalculation                    % calculate and store ExemplarIDI
zExemplarFrequencyCalculation

zWriteExemplarPDB(1)

zExemplarTable(1,0,0,1);
zExemplarTable(2,0,0,1);
zExemplarTable(3,0,0,1);
zExemplarTable(4,0,0,1);
zExemplarTable(5,0,0,1);
zExemplarTable(6,0,0,1);
zExemplarTable(7,0,0,1);
zExemplarTable(8,0,0,1);
zExemplarTable(9,0,0,1);
zExemplarTable(10,0,0,1);
zExemplarTable(11,0,0,1);
zExemplarTable(12,0,0,1);

zExemplarTable(13,0,0,1);
zExemplarTable(14,0,0,1);
zExemplarTable(15,0,0,1);

zExemplarTablesExcel

zExemplarTable(1,0,1,2);
zExemplarTable(2,0,1,2);
zExemplarTable(3,0,1,2);
zExemplarTable(4,0,1,2);
zExemplarTable(5,0,1,2);
zExemplarTable(6,0,1,2);
zExemplarTable(7,0,1,2);
zExemplarTable(8,0,1,2);
zExemplarTable(9,0,1,2);
zExemplarTable(10,0,1,2);
zExemplarTable(11,0,1,2);
zExemplarTable(12,0,1,2);

