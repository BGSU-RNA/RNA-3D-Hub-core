% xAnnotateWithKnownMotifs(File) finds all instances of motifs from the motif
% library in the folder MotifLibrary in the given file(s)
 
function [File] = xAnnotateWithKnownMotifs(File,Verbose,WriteHTML,TheseMotifs)

if nargin < 2,
  Verbose = 0;
end

if nargin < 3,
  WriteHTML = 0;
end

if nargin < 4,
  Motif = dir([pwd filesep 'MotifLibrary' filesep 'LIB*.mat']);
else
  for i = 1:length(TheseMotifs),
    Motif(i).name = TheseMotifs{i};
  end
end

if strcmp(class(File),'char'),
  Filename = File;
  File = zAddNTData(Filename,0,[],Verbose);
end

UsingLibrary = 1;                                % so FR3D just searches

for ff = 1:length(File),
  FN        = upper(File(ff).Filename);
  Filenames = {File(ff).Filename};
  if isfield(File(ff),'Motifs'),
    MotNum = length(File(ff).Motifs);              % original number of motifs
  else
    MotNum = 0;
  end
  File(ff).Nucl(File(ff).NumNT+1).Motif = [];

  if WriteHTML > 0,

  LText{1} = ['<a href = "index.html">Return to FR3D home page for ' FN '</a><br>'];
  LText{2} = ['<a href = "' FN '_interactions.html">List of all pairwise interactions in ' FN '</a><br>'];
  LText{3} = ['<a href = "' FN '_basepairs.html">List of basepair interactions in ' FN '</a><br>'];
  LText{4} = ['<a href = "' FN '_stacking.html">List of stacking interactions in ' FN '</a><br>'];
%  LText{4} = ['<a href = "' FN '_basephosphate.html">List of base-phosphate interactions in ' FN '</a><br>'];
  LText{5} = ['<a href = "' FN '_motifs.html">List of motifs found in ' FN '</a><br>'];
  LText{6} = ['<a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=' FN '">PDB entry for ' FN '</a><br>'];
  LText{7} = ['<a href="../">Return to list of analyzed structures</a><br>'];
  LText{8} = ['<a href="../../basepairs">Basepair catalog</a><br>'];
  LText{9} = ['<a href="../../MotifLibrary/index.html">FR3D motif library</a><br>'];
  LText{10} = ['<a href="../../index.html">FR3D home page</a><br>'];
  LText{11} = ['<a href="http://rna.bgsu.edu">BGSU RNA group home page</a><br><br>'];

  DN = [pwd filesep 'Web' filesep 'AnalyzedStructures' filesep FN];

  if ~(exist(DN) == 7),        % if directory doesn't yet exist
    mkdir(DN);
  end

  ffid = fopen([DN filesep FN '_motifs.html'],'w');
  fprintf(ffid,'<html>\n<title>%s motif list\n</title>\n', FN);
  fprintf(ffid,'<body>\n');
  fprintf(ffid,'<h1>Motif list for %s</h1>\n', FN);

  for L = 1:length(LText),
    if L ~= 5,
      fprintf(ffid,'%s\n', LText{L});
    end
  end

%  fid = fopen([pwd filesep 'Web' filesep 'AnalyzedStructures' filesep 'All' filesep FN '_motifs.txt'],'w');

%  fprintf(fid,'PDB_ID\tMotif_name\tNucleotide_1_PDB_Number\tNucleotide_1_Chain\tNucleotide_1_Sequence_Position\n');

  end

  for m = 1:length(Motif),

   if ((Motif(m).name(9) == '_') && (Motif(m).name(4) ~= 'U')) || (nargin == 4),

    if nargin < 4,
      load(['MotifLibrary' filesep Motif(m).name]);
      MotifNumber = Motif(m).name(1:8);
      MotifName   = Motif(m).name(10:end);
    else
      load(['SearchSaveFiles' filesep Motif(m).name]);
      MotifNumber = 'x';
      MotifName   = Motif(m).name;
    end 

    if Verbose > 0,
      fprintf('\nSearching %s for %s\n', FN, Motif(m).name);
    end

    Query = Search.Query;

    clear Search
    clear Candidates

    xFR3DSearch

    [s,t] = size(Candidates);

    % the following line assumes that the GU packing motif has the text
    % 'GU_packing' in its name and that the GU pair are nucleotides 1 and 2

    if ~isempty(strfind(Motif(m).name,'GU_packing')), % GU packing motif
      t = 3;                                          % skip last nucleotide
    end

    File(ff).Motifs(m+MotNum).Name  = Motif(m).name;
    File(ff).Motifs(m+MotNum).Count = s;

    if s > 0,                                    % at least one candidate
      for c = 1:s,                               % go through candidates

       for n = 1:(t-1),                          % go through indices of cands
        clear Mot
        Mot.Name = strrep(MotifName,'.mat','');
        Mot.Index = m+MotNum;
        Mot.Number = MotifNumber;
        Mot.Indices = Candidates(c,1:(t-1));
        if isempty(File(ff).Nucl(Candidates(c,n)).Motif),
          L = 0;
        else
          L = length(File(ff).Nucl(Candidates(c,n)).Motif);
        end
        File(ff).Nucl(Candidates(c,n)).Motif(L+1) = Mot;
       end
      end

      [y,i] = sort(Candidates(:,1));               % sort by 1st nucleotide #
      Search.Candidates = Candidates(i,:);         % re-order candidates
      Search.Discrepancy = Search.Discrepancy(i);

      if WriteHTML > 0,

      fprintf(ffid,'<a name=%s>\n',MotifNumber);
      fprintf(ffid,'<h2><a href="%s">%s</a></h2>\n',['http://rna.bgsu.edu/FR3D/MotifLibrary/' MotifNumber '/instances.html'], strrep(Motif(m).name,'.mat',''));

      fprintf(ffid,'<pre>\n');
      Text = xListCandidates(Search,Inf,6,{MotifNumber});
      for c = 4:length(Text),
        fprintf(ffid,'%s\n',Text{c});
      end
      fprintf(ffid,'</pre>\n');

      end
    end
   end
  end

  if WriteHTML > 0,
    fclose(ffid);
  end
end


return

ff = 1;
for n = 1:length(File(ff).NT),
  if ~isempty(File(ff).Nucl(n).Motif),
    for m = 1:length(File(ff).Nucl(n).Motif),
      fprintf('%4d %s %2d %s\n', n, File(ff).Nucl(n).Motif(m).Number, File(ff).Nucl(n).Motif(m).Index, File(ff).Nucl(n).Motif(m).Name);
    end
  end
end
