
% text after the percent symbol are comments, ignored by Matlab

% Filenames = {'1s72','2avy'};      % files to be searched
Filenames = 'Nonredundant_2008_02_21_list';  % files to be searched
%Filenames = {'1s72','2avy'};      % files to be searched

GUIactive = 1;                 % so xFR3DSearch will use this Query
Letters = 'ACGU';                 % makes it easy to convert numbers to letters
Verbose = 1;                      % tell xFR3DSearch to give information about progress

Edges = {'tWH','tSS','cHW','tHS'};   % possible interactions between bases 
                                   % Edges is a "cell array"

Edges = {'cWW','cHW'};   % possible interactions between bases 

%for dd = 1:length(Edges),       % first pairwise interaction
for dd = 1:1,       % first pairwise interaction
 for ee = 2:length(Edges),      % second pairwise interaction
  for ii = 1:1,                      % code of first letter in the triple
   for jj = 3:3,                     % code of second letter in the triple
    for kk = 3:3,                    % code of third letter in the triple
     clear Query                   % remove any previous value of this variable
     clear Search
     
     % --------------------- set up the criteria for a symbolic search

     Query.Mask           = [Letters(ii) Letters(jj) Letters(kk)];
%Query.Mask = 'NNN';
     Query.Edges{1,2}     = Edges{dd};
     Query.Edges{2,3}     = Edges{ee};
     Query.Geometric      = 0;               % symbolic search
     Query.Description    = [Query.Mask ' ' Edges{dd} ' ' Edges{ee} ' triple symbolic'];
     Query.Name           = Query.Description;
     
     Query = xConstructQuery(Query);         % preliminary calculations

     xFR3DSearch                          % conduct a symbolic search

     if isfield(Search,'Candidates'),
       [L,t] = size(Search.Candidates);     % L is number of candidates
     else
        L = 0;
     end
     
     if L > 0,

       % --------------------- set up the criteria for a geometric search

       xListCandidates(Search);
       save(['SearchSaveFiles' filesep Search.SaveName], 'Search');
       
       clear Query
       Query.Mask           = [Letters(ii) Letters(jj) Letters(kk)];
       Query.Geometric      = 1;
       f = Search.Candidates(1,4);              % file number of 1st candidate

       for m = 1:3,
         Query.NTList{m} = [Search.File(f).NT(Search.Candidates(1,m)).Number '_' Search.File(f).NT(Search.Candidates(1,m)).Chain];
       end
       Query.Filename       = Search.File(f).Filename;

       Query.Edges{1,2} = 'cp';
       Query.Edges{2,3} = 'cp';

       Query.Edges{1,2} = 'cp ncp';
       Query.Edges{2,3} = 'cp ncp';

       Query.DiscCutoff     = 0.5;              % maximum discrepancy
       Query.Description    = [Query.Mask ' ' Edges{dd} ' ' Edges{ee} ' triple geometric'];
       Query.Name           = Query.Description;
       [File,QIndex] = zAddNTData(Query.Filename,0,File);  
                                             % load data for Query, if needed
       Query = xConstructQuery(Query,File(QIndex)); % preliminary calculations

       xFR3DSearch                       % conduct geometric search

       xListCandidates(Search);
       save(['SearchSaveFiles' filesep Search.SaveName], 'Search');
       xDisplayCandidates(File,Search);
     end
    end
   end
  end
 end
end



