% zAnalyzedFilesHTML produces an .html file listing basepairing and base
% stacking interactions for each molecule in File

% To run from the command line, matlab -r zAnalyzedFilesHTML('1s72',path)
% At least I think that will work
% It infers the extension .pdb

function [void] = zAnalyzedFilesHTML(File,datapath)

path(path,[pwd filesep 'FR3DSource']);

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char') || strcmp(class(File),'cell'),
  Filename = File;
  File = zAddNTData(Filename,0,[],1);
end

zBackboneCodes;                           % load conformation codes

for f = 1:length(File),
 if length(File(f).NT) > 1,
  clear chainoffset

  FN = upper(File(f).Filename);
  PDBFN = File(f).PDBFilename;

  Vers = num2str(File(f).ClassVersion);

  DataHeader1 = sprintf('# PDB_ID_FR3D_Version_%s\tSuitename_version_suitename.0.3.070628\tDangle_version_chiropraxis.0.63.070705',Vers);

  DataHeader2 = sprintf('PDB_ID\tInteraction\tNucleotide_1_Base\tNucleotide_1_PDB_Number\tNucleotide_1_Chain\tNucleotide_1_Sequence_Position\tNucleotide_2_Base\tNucleotide_2_PDB_Number\tNucleotide_2_Chain\tNucleotide_2_Sequence_Position');

  LText{1} = ['<a href = "index.html">Return to FR3D home page for ' FN '</a><br>'];
  LText{2} = ['<a href = "' FN '_interactions.html">List of all pairwise interactions in ' FN '</a><br>'];
  LText{3} = ['<a href = "' FN '_basepairs.html">List of basepair interactions in ' FN '</a><br>'];
  LText{4} = ['<a href = "' FN '_stacking.html">List of stacking interactions in ' FN '</a><br>'];
  LText{5} = ['<a href = "' FN '_base_phosphate.html">List of base-phosphate and base-ribose interactions in ' FN '</a><br>'];
  LText{6} = ['<a href = "' FN '_backbone_connectivity.html">List of backbone connectivity relations in ' FN '</a><br>'];
  LText{7} = ['<a href = "' FN '_backbone_conformation.html">List of backbone conformations found in ' FN '</a><br>'];
  LText{8} = ['<a href = "' FN '_motifs.html">List of motifs found in ' FN '</a><br>'];
  LText{9} = ['<a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=' FN '">PDB entry for ' FN '</a><br>'];
  LText{10} = ['<a href="../">Return to list of analyzed structures</a><br>'];
  LText{11} = ['<a href="../../basepairs">FR3D Basepair catalog</a><br>'];
  LText{12} = ['<a href="../../BasePhosphates">FR3D Base-phosphate interactions</a><br>'];
  LText{13} = ['<a href="../../MotifLibrary/index.html">FR3D motif library</a><br>'];
  LText{14} = ['<a href="../../index.html">FR3D home page</a><br>'];
  LText{15} = ['<a href="http://rna.bgsu.edu">BGSU RNA group home page</a><br><br>'];

  HText{1} = ['FR3D classification version ' num2str(File(f).ClassVersion) ' ' date];

  HText{2} = '<p>Basepairing follows the paper Leontis, Stombaugh, Westhof Nucleic Acids Research <b>30</b> No. 16, August 15, 2002.  See the <a href="../../basepairs">FR3D Basepair catalog</a>.  Basepairs are either <i>cis</i> or <i>trans</i>, denoted c and t below.  Each base can use one of three edges, the Waston-Crick edge (W), the Hoogsteen edge (H) or the Sugar edge (S).  Basepairs listed below indicate which base is using which edge.  For example, a line reading A108 G130 tSH would mean that A108 is using its Sugar Edge and G130 is using its Hoogsteen edge.  In the case that both bases use the same edge, a capital letter indicates which base uses the edge in the dominant way.  For perfectly symmetric basepairs such as AU cWW, the capital and lowercase letters are irrelevant.  Bifurcated basepairs are indicated by the text bif.  It does not indicate which base uses which edge.';

  HText{3} = '<p>Base stacking is divided into three subcategories, according to the faces used by each base.  The faces are named this way:  Imagine a base making a Watson-Crick basepair in a helical context.  The side of the base that faces toward the 3'' end of the strand is called the 3 face, while the side that faces the 5'' end is called the 5 face.  The stacking one finds in a helix is thus refered to as s35, meaning that the first base uses the 3 face and the second uses the 5 face.  Two other types of stacking (s33 and s55) occur in other contexts.';

  HText{4} = '<p>Each pair is listed twice.  For instance, if A108 G130 tSH is listed, then so is G130 A108 tHS.  The order in which the edges are listed still corresponds to the base which is using that edge.  Similarly with stacking.  The chain is indicated in parentheses.';

  HText{5} = '<p>Starting from hairpin loops, cWW interactions are classified as being nested if they do not cross any previously established nested cWW pairs.  The last column indicates, for all interactions, the number of nested cWW pairs crossed by the interaction (which may be visualized on the circular diagram of pairwise interactions).';

  HText{6} = '<p>Classification of basepairs and base stacking is done using the program FR3D developed by the <A href="http://rna.bgsu.edu">Bowling Green State University RNA group</a>, including Neocles Leontis, Craig L. Zirbel, Jesse Stombaugh, Ali Mokdad, Michael Sarver, and Anton Petrov.';

  HText{7} = '<p>Base-phosphate interactions are described in <a href="http://nar.oxfordjournals.org/cgi/reprint/gkp468v1?ijkey=Gs0yMo4zWz3Po7X&keytype=ref">this paper</a> by the BGSU RNA group.  See the catalog of <a href="../../BasePhosphates">FR3D Base-phosphate interactions</a>.  Base-ribose interactions are annotated using the same method.  Backbone connectivity relations are denoted c35 or c53.  c35 indicates that the first nucleotide uses its O3'' atom to connect to the phosphorus of the second nucleotide and then to the O5'' atom of the second nucleotide.  c53 is used when the order of the nucleotides is reversed.';

  HText{8} = '<p>Backbone conformations are two-character codes which apply to the portion of the backbone shared by the two listed nucleotides.  They are only listed once, with the nucleotides in increasing order.  Nucleotides missing a base in the PDB file are omitted by FR3D.  Conformations are calculated by Dangle and Suitename from the <A href="http://kinemage.biochem.duke.edu/index.php">Richardson lab</A> at Duke University.'; 

  HText{9} = '<p>Please write to Craig L. Zirbel at <img src="http://www-math.bgsu.edu/z/cz.gif" align="absbottom"> with questions or comments.<hr>';

  HText{10} = '<pre>';

  % -------------------------------------------------- Set paths

  warning off

  if nargin > 1,

    mypath = [datapath filesep FN filesep];
    datapath = [datapath filesep 'All' filesep];

  else

    mypath = [pwd filesep 'Web' filesep 'AnalyzedStructures'];

    mkdir(mypath, FN);
    warning on

    mypath = [mypath filesep FN filesep];

    datapath = [pwd filesep 'Web' filesep 'AnalyzedStructures' filesep 'All' filesep FN filesep];

    mkdir(datapath)

  end


  % ------------------------------------------- Write chains and sequences

%  fid = fopen([datapath FN '_chain_sequence_FR3D.txt'],'w'); % open for writing
%  fprintf(fid,'%s\n',DataHeader1);
%  fclose(fid);

  Chain = cat(2,File(f).NT.Chain);
  U     = unique(Chain);
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    chainoffset(u) = min(i);                    % first index in this chain
  end

  % ------------------------------------------- Produce interaction list

  c = 1;                                    % counter for interactions

  IText{1} = '';
  DText{1} = '';
  InterType = [];

  E   = File(f).Edge;
  BPh = File(f).BasePhosphate;
  BR  = File(f).BaseRibose;
  BC  = File(f).Covalent;                    % covalent connections

  for i = 1:File(f).NumNT,
    N1 = File(f).NT(i);

    % ------------------------------------- Find basepairing, stacking
    j = find(E(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', full(File(f).Crossing(i,j(k))));
      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %7s - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zEdgeText(File(f).Edge(i,j(k)),0,N1.Code,N2.Code),r);

      u  = find(U==N1.Chain);              % which chain i is in
      ii = i - chainoffset(u) + 1;        % position of i in chain u
      u  = find(U==N2.Chain);              % which chain j(k) is in
      jj = j(k) - chainoffset(u) + 1;     % position of j(k) in chain u

      T = zEdgeText(File(f).Edge(i,j(k)),0,N1.Code,N2.Code);

      DText{c} = sprintf('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n', File(f).Filename, T, N1.Base, N1.Number, N1.Chain, ii, N2.Base, N2.Number, N2.Chain, jj);

      InterType(c) = abs(File(f).Edge(i,j(k)));

      c = c + 1;
    end

    % ------------------------------------- Find base phosphate interactions
    j = find(BPh(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', full(File(f).Crossing(i,j(k))));

      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %7s - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zBasePhosphateText(BPh(i,j(k)),1), r);

      u  = find(U==N1.Chain);              % which chain i is in
      ii = i - chainoffset(u) + 1;        % position of i in chain u
      u  = find(U==N2.Chain);              % which chain j(k) is in
      jj = j(k) - chainoffset(u) + 1;     % position of j(k) in chain u

      T = zBasePhosphateText(BPh(i,j(k)),1);

      DText{c} = sprintf('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n', File(f).Filename, T, N1.Base, N1.Number, N1.Chain, ii, N2.Base, N2.Number, N2.Chain, jj);

      if abs(BPh(i,j(k))) > 100,
        InterType(c) = abs(BPh(i,j(k)));      % near interaction
      elseif i == j(k),
        InterType(c) = 200.1;                 % self interaction
      else
        InterType(c) = 200;                   % non-self interaction
      end

      c = c + 1;
    end

    % ------------------------------------- add base phosphate interactions
    j = find(BR(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', full(File(f).Crossing(i,j(k))));

      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %7s - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zBaseRiboseText(BR(i,j(k)),1), r);

      u  = find(U==N1.Chain);              % which chain i is in
      ii = i - chainoffset(u) + 1;        % position of i in chain u
      u  = find(U==N2.Chain);              % which chain j(k) is in
      jj = j(k) - chainoffset(u) + 1;     % position of j(k) in chain u

      T = zBaseRiboseText(BR(i,j(k)),1);

      DText{c} = sprintf('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n', File(f).Filename, T, N1.Base, N1.Number, N1.Chain, ii, N2.Base, N2.Number, N2.Chain, jj);

      if abs(BR(i,j(k))) > 100,
        InterType(c) = abs(BR(i,j(k)));      % near interaction
      elseif i == j(k),
        InterType(c) = 200.1;                 % self interaction
      else
        InterType(c) = 200;                   % non-self interaction
      end

      c = c + 1;
    end

    % ---------------------------------- Find backbone connectivity relations
    j = find(BC(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', 0);

      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %6s  - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zBackboneContinuityText(BC(i,j(k))),r);

      u  = find(U==N1.Chain);              % which chain i is in
      ii = i - chainoffset(u) + 1;        % position of i in chain u
      u  = find(U==N2.Chain);              % which chain j(k) is in
      jj = j(k) - chainoffset(u) + 1;     % position of j(k) in chain u

      T = zBackboneContinuityText(BC(i,j(k)));

      DText{c} = sprintf('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n', File(f).Filename, T, N1.Base, N1.Number, N1.Chain, ii, N2.Base, N2.Number, N2.Chain, jj);

      InterType(c) = 300;                 % code for later

      c = c + 1;
    end

    % ---------------------------------- Find backbone conformation
    j = find(File(f).Backbone(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', 0);

      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %6s  - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, Codes{File(f).Backbone(i,j(k))},r);

      u  = find(U==N1.Chain);              % which chain i is in
      ii = i - chainoffset(u) + 1;        % position of i in chain u
      u  = find(U==N2.Chain);              % which chain j(k) is in
      jj = j(k) - chainoffset(u) + 1;     % position of j(k) in chain u

      T = Codes{File(f).Backbone(i,j(k))};

      DText{c} = sprintf('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n', File(f).Filename, T, N1.Base, N1.Number, N1.Chain, ii, N2.Base, N2.Number, N2.Chain, jj);

      InterType(c) = 400;                 % code for later

      c = c + 1;
    end
  end

% -------------------------------------------------- Write index.html file

  fid = fopen([mypath 'index.html'],'w');

  

  fprintf(fid,'<html>\n<head>\n<title>%s analysis by FR3D</title>\n<script src="../../jmol/Jmol.js" type="text/javascript"></script></head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D analysis of ' PDBFN '</h1>']);

  for L = 2:length(LText),
    fprintf(fid,'%s\n', LText{L});
  end

  Diagram = ['<a href="' FN '_circular_diagram.pdf"> <img src="' FN '_circular_diagram.png" alt="Click for high resolution pdf" width = "615" > </a>'];

  DiagramText = '<p>Circular basepair diagram in Nussinov style.  Click the diagram for a high-resolution PDF version.';

%  Dark blue chords indicate nested Watson-Crick basepairs, cyan indicates nested non-Watson-Crick basepairs, red indicates non-nested Watson-Crick basepairs, green indicates non-nested non-Watson-Crick basepairs, and yellow indicates long-range stacking interactions.

  fprintf(fid,'<table border="1">\n<tr>\n<th>Circular interaction diagram</th>\n<th>RNA 3D structure</th>\n</tr>\n');

  fprintf(fid,'<tr>\n<td>\n%s\n</td>\n',Diagram);

  fprintf(fid,'<td>\n<script type=''text/javascript''>\n jmolInitialize(''../../jmol'');\n jmolApplet(600, ''load %s_RNA.pdb;spacefill off'',''3D structure'');\n </script>\n</td>\n</tr>\n', FN);

  fprintf(fid,'<tr>\n<td>\n%s\n</td>\n', DiagramText);

  fprintf(fid,'<td><a href="%s_RNA.pdb">Click here</a> to download the RNA coordinate file (first model only when multiple models exist).\n</td>\n</tr>\n</table>\n',FN);

  fprintf(fid,'<b>Resolution: </b>%7.1f<br>\n', File(f).Info.Resolution);
  fprintf(fid,'<b>Descriptor: </b>%s<br>\n', File(f).Info.Descriptor);
  fprintf(fid,'<b>Experimental technique: </b>%s<br>\n', File(f).Info.ExpTechnique);
  fprintf(fid,'<b>Release Date: </b>%s<br>\n', File(f).Info.ReleaseDate);
  fprintf(fid,'<b>Author: </b>%s<br>\n', File(f).Info.Author);
  if ~isempty(File(f).Info.Source),
    fprintf(fid,'<b>Biological source: </b>%s<br>\n', File(f).Info.Source);
  end

  fprintf(fid,'</html>\n');

  fclose(fid);

  % --------------------------------------------- Write FN_interactions file  

  fid = fopen([mypath FN '_interactions.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s interactions from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D list of pairwise interactions in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 2,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find((InterType < 30) + (InterType >= 200));  % exclude near pairs, BPh

  for i = 1:length(k),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

% datapath

  if ~exist([datapath FN],'dir'), mkdir([datapath FN]); end % Anton 6/3/2011

  fid = fopen([datapath FN '_interactions_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

  % ----------------------------------------- Write FN_near_interactions file  

  fid = fopen([mypath FN '_near_interactions.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s near interactions from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D list of near pairwise interactions in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 2,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find((InterType > 100) .* (InterType < 200) .* (fix(InterType) ~= 114));  % exclude near pairs, BPh

  for i = 1:length(k),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

  fid = fopen([datapath FN '_near_interactions_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

  % ----------------------------------------------- Write FN_basepairs file  

  fid = fopen([mypath FN '_basepairs.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s basepairs from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D list of basepair interactions in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 3,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find(InterType < 14);                    % exclude Rib pairs

  for i = 1:length(k),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

  fid = fopen([datapath FN '_basepairs_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

% ----------------------------------------------- Write FN_stacking file  

  fid = fopen([mypath FN '_stacking.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s stacking interactions from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D list of stacking interactions in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 4,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find((InterType > 19) .* (InterType < 100));

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

  fid = fopen([datapath FN '_stacking_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

% --------------------------------------------- Write FN_base_phosphate file  

  fid = fopen([mypath FN '_base_phosphate.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s base-phosphate and base-ribose interactions from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D list of base-phosphate and base-ribose interactions in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 5,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:(length(HText)-2),
    fprintf(fid,'%s\n',HText{i});
  end

  fprintf(fid,'<p>Non-self interactions listed first, if any, then self interactions.\n');

  for i = (length(HText)-1):length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find(InterType == 200);               % non-self interactions

  c = 1;

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',c,IText{k(i)});
    c = c + 1;
  end

  k = find(InterType == 200.1);               % self interactions

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',c,IText{k(i)});
    c = c + 1;
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

  fid = fopen([datapath FN '_base_phosphate_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  k = find(fix(InterType) == 200);               % all BPh and BR interactions
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

% --------------------------------------- Write FN_backbone_connectivity file  

  fid = fopen([mypath FN '_backbone_connectivity.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s backbone connectivity relations from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D backbone connectivity relations in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 6,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find(InterType == 300);

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

  fid = fopen([datapath FN '_backbone_connectivity_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

% --------------------------------------- Write FN_backbone_conformation file  

  fid = fopen([mypath FN '_backbone_conformation.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s backbone conformation relations from Richardson lab programs</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>Backbone conformation relations in ' FN '</h1>']);
  fprintf(fid,'As computed by Dangle and Suitename from the <A href="http://kinemage.biochem.duke.edu/index.php">Richardson lab</A> at Duke University</p>\n');

  for L = 1:length(LText),
    if L ~= 7,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find(InterType == 400);

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

  fid = fopen([datapath FN '_backbone_conformation_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);


  clear IText HText DText
 end
end


