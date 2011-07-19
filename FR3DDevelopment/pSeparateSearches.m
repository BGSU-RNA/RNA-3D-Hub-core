% pSeparateSearches.m takes the searches in MotifLibrary and writes them out as individual .fasta files.  We no longer need individual .mat files, but that code is at the bottom.

% ------------------------------------------- Read names of model files

disp('Make sure the Matlab current folder has a MotifLibrary in it');

if ~exist('loopType'),
  loopType = 'HL';
  %loopType = 'IL';
end

% ------------------------------------------- Keep only loopType files

keep = [];
Filenames = dir(['MotifLibrary' filesep]);
for m = 1:length(Filenames),
    if (length(Filenames(m).name) > 2),
      if (Filenames(m).name(1:2) == loopType), 
          keep(m) = 1;
          Filenames(m).name = strrep(Filenames(m).name,'.mat','');
      end
    end 
end
Filenames = Filenames(find(keep));           % keep only these models

% -------------------------------------------- Loop through models, write fasta

nameid = fopen(['Sequences' filesep loopType '_SeparatedSequences.txt'],'w');

for i = 1:length(Filenames),                   % go through all models
  MN = Filenames(i).name;
  FN = ['MotifLibrary' filesep MN '.mat'];     % .mat file
  fullSearch = load(FN,'Search','-mat');       % load .mat file
  [m,n] = size(fullSearch.Search.Candidates);  % m candidates, n nucleotides
  n = n - 1;
  Text = xFASTACandidates(fullSearch.Search.File,fullSearch.Search,0,MN);

  for j = 1:m,                                   % loop through candidates
    f = fullSearch.Search.Candidates(j,n+1);     % file number of candidate
    i = fullSearch.Search.Candidates(j,1:n);     % indices of nucleotides
    E = fullSearch.Search.File(f).Edge(i,i);     % interactions between them

    Sig = zMotifSignature(E);                    % motif signature

    if m >= 100,
      FN = [MN '_' sprintf('%3d',j) '_' Sig '.fasta'];
    elseif m >= 10,
      FN = [MN '_' sprintf('%2d',j) '_' Sig '.fasta'];
    else
      FN = [MN '_' sprintf('%1d',j) '_' Sig '.fasta'];
    end
    FN = strrep(FN,' ','0');

    fprintf(nameid,'%s\n', FN);                  % write filename to list

    fid = fopen(['Sequences' filesep FN],'w');
    fprintf(fid,'%s\n',Text{2*j-1});
    fprintf(fid,'%s\n',Text{2*j});
    fclose(fid);

  end
end

fclose(nameid);

break


if 0 > 1,
for i = 1:length(Filenames),                   % go through all models
  MN = Filenames(i).name;
  FN = ['MotifLibrary' filesep MN '.mat'];     % .mat file
  fullSearch = load(FN,'Search','-mat');       % load .mat file
  [m,n] = size(fullSearch.Search.Candidates);  % m candidates, n-1 nucleotides
  for j = 1:m,                                 % loop through candidates
      f = fullSearch.Search.Candidates(j,n);           % file numbers of motifs
      Search.Candidates = fullSearch.Search.Candidates(j,:);
      Search.Candidates(n) = 1;                % file number 1, only 1 instance
      Search.File = fullSearch.Search.File(f); % File it came from
      Search.Discrepancy = fullSearch.Search.Discrepancy(j);
      Search.Query = fullSearch.Search.Query;
      Search.origmatfilename = fullSearch.Search.origmatfilename;
      Search.Signature = fullSearch.Search.Signature;
      Search.matfilename = fullSearch.Search.matfilename;
      Search.Truncate = fullSearch.Search.Truncate; % might be for 1st instance?
      Search.ownsequencefasta = fullSearch.Search.ownsequencefasta;
      Search.modelfilename = fullSearch.Search.modelfilename;
      oFN = ['MotifLibrary' filesep 'Separated' filesep MN '_' int2str(j) '.mat'];
      save(oFN,'Search','-mat')
  end
end
end
% now cd to Motifs\Separated and run pMakeModelsFromLibrary
