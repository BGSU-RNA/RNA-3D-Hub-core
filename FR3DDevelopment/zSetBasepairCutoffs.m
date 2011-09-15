% zSetBasepairCutoffs finds true and near instances of basepairs and provides diagnostics to help modify cutoffs


%%% load NRTemp


NRList = 'Nonredundant_4A_2010-05-19_list';
NRList = 'Nonredundant_4A_2011-01-07_list';
NRList = 'Nonredundant_4A_2011-07-23_list';

Verbose = 1;               % Verbose = 1 tells it to show distance graphs
LMax    = 500;             % maximum number of pairs to consider in each class

% Pair codes:  1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU

pcodes = [9 11 13 14 15 16 1 5 6 7];    % pair codes to work on

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

if ~exist('Search'),
  if ~exist('File'),                           % if no molecule data is loaded,
    [File,SIndex] = zAddNTData(NRList,0,[],1);   % load PDB data
    File = File(SIndex);
  else
    [File,SIndex] = zAddNTData(NRList,0,File,1); % add PDB data if needed
    File = File(SIndex);
  end                       
  Search.File = File;
  clear File
end

if isfield(Search.File,'AA'),
  Search.File = rmfield(Search.File,'AA');      % removed 700 MB in July 2011!
end

Search.SaveName = 'BasepairTemp';

% loop through paircodes and computer classifications ----------------------

for j = 1:length(pcodes),            % run through all pair codes specified
 pc = pcodes(j);                     % current paircode

 CLE = CL(:,1,pc);                   % class limits for this paircode
 CLE = CLE(find(CLE));               % leave out empty entries
 CLE = CLE(find(abs(CLE) < 20));     % leave out stacking

 for row = 1:length(CLE),            % run through classes for this paircode

   % specify criteria for selection of pairs ----------------------------------
   % select pairs using selection criteria -----------------------------------

   if abs(CLE(row)) < 20,                  % if a basepair,
     Decimal = 1;                          % treat subcategories separately
   else
     Decimal = 0;                          % group subcategories together
   end

   Decimal = 0;



   Param = [];
   Param(1,1) = CLE(row);
   Param(1,2) = pc;

   [ListTrue,Class] = zFindPairs(Search.File,Param,Decimal);   % search the non-redundant list

   Param(1,1) = sign(CLE(row)) * (100 + abs(CLE(row)));  % near pairs too!
   Param(1,2) = pc;

   ListNear = zFindPairs(Search.File,Param,Decimal);   % search the non-redundant list

   if length(ListTrue) > 0,
     if length(ListNear) > 0,
       List = [ListTrue zeros(size(ListTrue(:,1))); ListNear ones(size(ListNear(:,1)))];
       Class = [Class; zeros(length(ListNear(:,1)),1)];
     else
       List = ListTrue;
       Class = Class;
     end
   elseif length(ListNear) > 0,
     List = ListNear;
     Class = zeros(length(ListNear(:,1)),1);
   else
     List = [];
     Class = [];
   end

   if length(List) > 0,               % instances of this category are found

    fprintf('Found %5d true instances of %2s %4s class %5.0f\n', length(ListTrue), Pairs{pc}, zEdgeText(CLE(row),Decimal,pc), CLE(row));

    fprintf('Found %5d near instances of %2s %4s class %5.0f\n', length(ListNear), Pairs{pc}, zEdgeText(CLE(row),Decimal,pc), CLE(row));

    if length(List(:,1)) > LMax,
      fprintf('Randomly choosing %d of them\n', LMax);
      r = rand(size(List(:,1)));
      [y,i] = sort(r);                    % randomize the order of the instances
      SList = List(i(1:min(end,LMax)),:);  % keep first LMax of them only
    else
      SList = List;
    end

    [PD,ID,i,k] = zOrderPairs(Search.File,SList,LMax,Verbose); % order by centrality

    PD = PD(i,i);
    List = List(i,:);                   % order by discrepancy

    for i = 1:length(List(:,4)),        % use diagonal to tell near and subcategory
      if List(i,4) == 0,
        
        PD(i,i) = abs(abs(Class(i)) - fix(abs(Class(i))));
        if PD(i,i) > 0,
          PD(i,i) = PD(i,i) + 0.1;      % emphasize the color a bit more
        end
      else
        PD(i,i) = max(max(PD));         % max will appear blue
      end
    end

    Search.Candidates = List(:,1:3);
    Search.Disc = PD;
    Search.DiscComputed = ones(1,length(List(:,1)));
    Search.Query.Geometric = 1;
    Search.Query.NumNT = 2;
    Search.Query.DiscCutoff = 0.5;
    Search.Query.RelCutoff = 0.5;
    Search.Discrepancy = PD(:,1);

    if (length(List(:,1)) > 2) && (Verbose > 0),
      figure(8)
      clf
      q = zClusterGraph(PD);
      colormap('default');
      map = colormap;
      map = map((end-8):-1:8,:);
      colormap(map);
      caxis([0 1]);
      colorbar('location','eastoutside');
      hold on
      j = find(q == i(1));
      plot(j+0.5,j+0.5,'w*');
      m = find(q == k(1));
      plot(m+0.5,m+0.5,'wo');

      title(['Discrepancy between instances of ' Pairs{pc} ' ' zEdgeText(CLE(row),Decimal,pc)]);
      xlabel('Centroids marked (Geo discrep: white star, IDI: white circle)');
      ylabel('Blue on diagonal means non coplanar');
      drawnow
    end

    if 10 > 1,
      figure(1)
      close
      figure(1)
%      xScatterPairs(Search,1,2);
%      xDisplayCandidates(Search.File,List(:,1:3));
      xDisplayCandidates([],Search);
      rotate3d on
    end

    end
  end
end
