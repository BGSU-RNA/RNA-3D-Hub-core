% zUpdatePDBDatabase reads a file holding the current list of PDB files and metadata, reads all 3D structure files, shares information between these two sources, saves PDBInfo.mat, then finds non-redundant lists and writes them out, then finds exemplar basepairs. 

% zUpdatePDBDatabase('2010-11-12',600);
% zUpdatePDBDatabase('2010-11-12');
% zUpdatePDBDatabase('2010-12-08');
% zUpdatePDBDatabase('2010-12-15');
% zUpdatePDBDatabase('2011-01-07');

function [] = zUpdatePDBDatabase(reportdate,current,ReadCode)

    reportdate = strrep(reportdate,'.txt','');           % in case it was a file
    reportdate = strrep(reportdate,'report_','');        % with a long filename

    if nargin < 3,
      ReadCode = 0;
    end

    Verbose = 1;

    filename = ['report_' reportdate '.txt'];

    filename

    % load PDBInfo

    if nargin < 2,
      current = 1;
    end

    initial = current;

    if current == 1,
      fid = fopen(filename,'r');
      temp = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%n\t%s', 'delimiter', '\t');
      fclose(fid);

      t = cell(length(temp{1}),length(temp));
      n = zeros(length(temp{1}),3);
      for i = 1:length(temp{1})
        for j = 1:length(temp)
          if j == 7
            n(i,1) = temp{j}(i);         % store resolution
            t{i,j} = '';
          else
            t{i,j} = temp{j}{i};
          end
        end    
      end
    else
%       load PDBInfo
        load([pwd filesep 'FR3DSource' filesep 'PDBInfo.mat']);
    end

    fprintf('Updating %d PDB files using %s\n', length(t(:,1)), filename);

    % Columns should correspond to this list:
    % 1 A	Structure ID
    % 2 B	Descriptor / Structure title
    % 3 C	Experimental Technique / Exp. method
    % 4 D	Release Date
    % 5 E	Authors
    % 6 F	Keywords
    % 7 G	Resolution (Å)
    % 8 H	Source

    % We should add an additional column for NDB ID	
    % Later in this program, other columns are set in t and in n:
    % 11    Sequence of longest chain
        % 12    Best chains, as determined by zBestChains

        % In zFileRedundancy, column 10 is set
        % 10    PDB ID of structure which represents this one

    % -------------------------------- Create PDBInfo.mat

    t{1,12} = '';                     % make space to save base sequence
    n(1,5) = 0;                      % make space to save numeric values

    % --------------------------------- Remove duplicate entries - May 2010
    Keep = zeros(length(t(:,1)),1);             % which entries to keep

    Keep(1) = 1;
    first = 1;
    i = 2;

    while i <= length(t(:,1)),

      while (i <= length(t(:,1))) && strcmp(lower(t{first,1}),lower(t{i,1})),
                                                % first, i are duplicate entries
        if ~strcmp(t{first,2},t{i,2}),
          fprintf('Different descriptors for %s\n', t{first,1});
          fprintf('%s\n', t{first,2});
          fprintf('%s\n', t{i,2});
        end

        i = i + 1;

      end

    %  fprintf('Kept %s\n', t{i,1});

      Keep(i) = 1;
      first = i;
      i = i + 1;

    end

    j = find(Keep);
    t = t(j,:);
    n = n(j,:);

%     save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7
    save([pwd filesep 'FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7
    % ------------------------------------- Read each PDB file now, save data


    for i = current:length(t(:,1)),

      fprintf('Updating structure %d, which is %s ===============================\n', i, t{i,1});

      current = i;

      try
        File = zAddNTData(t{i,1},ReadCode,[],1);          % load file
      catch
        delete(['PrecomputedData' filesep t{i,1} '.mat']);
        File = zAddNTData(t{i,1},ReadCode,[],1);          % load file
      end

      if 0 > 1,
        File = zGetPDBInfo(File,n,t);
        File.Info
        zSaveNTData(File);
      end

      if 0 > 1,
        if length(File.NT) > 1,                    % if it has nucleotides,
          c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers
          File.Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
          File.Distance(1,1) = 10;                 % make sure one is non-zero
          d = sort(nonzeros(File.Distance));

          if d(min(10,length(d))) < 1,
            fprintf('%s might have overlapping nucleotides, loading it again\n',File.Filename);
            File = zAddNTData(t{i,1},4,[],1);          % load file
          end        
        end
      end

      if ~isempty(File.NT),                      % if it has nucleotides,

        E  = abs(triu(File.Edge));
        n(i,3) = nnz(zSparseRange(E,0,13));      % number of regular pairs
        n(i,2) = length(File.NT);                % store the number of NT

        LC = File.LongestChain;

        t{i,11} = cat(2,File.NT(LC(1):LC(2)).Base);    % bases in longest chain
        t{i,12} = File.BestChains;        % characters of the best chain(s)

        n(i,4) = length(t{i,11});         % number of nucleotides in longest chain
        n(i,5) = LC(1);                   % starting index of longest chain
        n(i,6) = LC(2);                   % end index of longest chain

        if Verbose > 1,
          fprintf('All      %s\n',cat(2,File.NT.Base));
          fprintf('Longest  %s\n',t{i,11});
        end
      end

      if mod(i,50) == 0 || i == length(t(:,1)),
%         save(['PDBInfo.mat'],'n','t'); % Matlab version 7
            save([pwd filesep 'FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7
      end

    end

%     save(['PDBInfo.mat'],'n','t'); % Matlab version 7
    save([pwd filesep 'FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

    if initial == 1,                     % we've gone through the whole list
      if exist('Dropboxroot'),
        save([DropboxRoot filesep 'FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7
      end
    end

    % ------------------------------------- Find a non-redundant list of PDBs

    zFileRedundancy_2(reportdate,t,n)

    zWriteHTMLFileList(reportdate)              % write out NR lists as HTML

    % ------------------------------------- Find new exemplars using NR list

    % zFindExemplars






    % OLD CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	 download a new report, save it as PDB_File_Information current-date.xls
    %	 and modify the following line:
    %	
    %	Date = '2010-05-19';
    %	
    %	PDBInfoName = ['PDB_File_Information ' Date '.xls'];

    %	[n,t] = xlsread(PDBInfoName);  % read Excel report on PDB files
    %	
    %	t = t(2:end,:);                  % remove the header row
