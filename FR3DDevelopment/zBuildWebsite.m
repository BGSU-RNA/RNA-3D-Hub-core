
WhatToWrite = [0 1 0 0 0];           % which pieces to write right now
WhatToWrite = [1 0 0 0 0];           % which pieces to write right now
WhatToWrite = [0 2 0 0 0];           % which pieces to write right now
WhatToWrite = [1 1 1 0 0];           % which pieces to write right now
WhatToWrite = [1 1 1 0 1];           % which pieces to write right now

% WhatToWrite(1)  - html pages
% WhatToWrite(2)  - circular diagram
% WhatToWrite(3)  - pdb file
% WhatToWrite(4)  - annotate with known motifs
% WhatToWrite(5)  - only work with new files, where no web directory exists

load PDBInfo

Names = {'2AW4','2J01','1S72','2AVY','1J5E','3BWP','1ZZN'};
Names = {'1HR2'};
Names = {'3BWP'};
%Names = zReadPDBListim('Nonredundant_2008_02_21_list');
Names = t(:,1);                        % names of files from PDB/NDB

WTW = WhatToWrite;

current = length(Names);
current = 1;

% current = find(ismember(t(:,1),'3PYO'));

for f = current:length(Names),  
%for f = length(Names):-1:1,  

  current = f;
  tim = cputime;

  if exist(['Web' filesep 'AnalyzedStructures' filesep Names{f}]) == 7,
    New = 0;
  else
    New = 1;
  end

  if New > 0 || WTW(5) == 0,

    fprintf('Writing HTML files for %s, file %d of %d\n', Names{f}, f, length(Names));

    ND = [pwd filesep 'Web' filesep 'AnalyzedStructures' filesep Names{f}];

    if New == 1,
      mkdir(ND);
    end

    File = zAddNTData(Names{f},0,[],1);              % load RNA data

    tim(end+1) = cputime;

    if isempty(File.NT),

      fprintf('No nucleotides found in this structure\n');

      fid = fopen([ND filesep 'index.html'],'w');
      fprintf(fid,'<html>\n');
      fprintf(fid,'No complete RNA nucleotides were found in this structure\n');
      fprintf(fid,'</html>\n');
      fclose(fid);

    else

      mypathroot = [pwd filesep 'Web' filesep 'AnalyzedStructures'];
      mypath     = [mypathroot filesep File.Filename filesep];

      % ----------- Make html to compare equivalent structures

        maxequiv = 21;
        k = find(ismember(t(:,10),t{f,10}));    % others in the class

        [y,i] = sortrows([n(k,3) ./ n(k,2) n(k,1)], [-1 2]);

        k = setdiff(k,f);
        k = [f; k];                            % put representative first

        fid = fopen([ND filesep t{f,1} '_compare_circular.html'],'w');
        fprintf(fid,'<html>\n');
        fprintf(fid,'<title>\n');
        fprintf(fid,'%s and %d equivalent structures\n',t{f,1},min(maxequiv,length(k))-1);
        fprintf(fid,'</title>\n');

        fprintf(fid,'<style>img{width:30%%;height:auto}</style>\n');

        fprintf(fid,'<body>\n');

        fprintf(fid,'Circular diagrams of %s and %d equivalent structures, in order of decreasing basepairs per nucleotide\n',t{f,1},min(maxequiv,length(k))-1);
        fprintf(fid,'<nobr>\n');
        for kk = 1:min(maxequiv,length(k)),
          fprintf(fid,'<a href="../%s/%s_circular_diagram.pdf"><img src="../%s/%s_circular_diagram.png" alt="Try another browser" width="2751" height="3181"></a>',t{k(kk),1},t{k(kk),1},t{k(kk),1},t{k(kk),1});
%          fprintf(fid,'<img src="../%s/%s_circular_diagram.png" alt="Click for high resolution pdf" width="2751" height="3181">',t{k(kk),1},t{k(kk),1});
        end
        fprintf(fid,'</nobr>\n');
        fprintf(fid,'</body>\n');
        fprintf(fid,'</html>\n');
        fclose(fid);

      % ----------------------------------------------- Write HTML files

      if WhatToWrite(1) > 0,
        zAnalyzedFilesHTML(File,mypathroot);             % make HTML
      end

      tim(end+1) = cputime;

      % ----------------------------------------------- Write PDB file

      if WhatToWrite(3) > 0,
        zWritePDB(File,[mypath File.Filename '_RNA.pdb']);
      end
      tim(end+1) = cputime;

      % ------------------------------------------- Annotate with known motifs

      if WhatToWrite(4) > 0,
        xAnnotateWithKnownMotifs(File,1,1);           % find and list motifs
      end
      tim(end+1) = cputime;

      % --------------------------------------------- Create circular diagram

      if WhatToWrite(2) > 0,
        clf
        zCircularDiagram(File,0.1,[1 1 1 1 1 1 1 1 1 1 1]);
        saveas(gcf,[mypath File.Filename '_circular_diagram.pdf'],'pdf');

        clf
        zCircularDiagram(File,1,[1 1 1 1 1 1 1 1 1 0 0]);

        tim(end+1) = cputime;

        FN = File.Filename;
        clear File

        print(gcf,'-dpng','-r600',[mypath FN '_circular_diagram.png']);
%        saveas(gcf,[mypath FN '_circular_diagram.png'],'png');

        [X,map] = imread([mypath FN '_circular_diagram.png']);
%        X = X(30:830,210:1030,:);
        X = X(120:3300,1250:4000,:);                % crop image
        imwrite(X,[mypath FN '_circular_diagram.png']);
        clear X
      end
      tim(end+1) = cputime;

      fprintf('Time taken:');
      fprintf(' %6.2f', diff(tim));
      fprintf(' Total: %6.2f', tim(end)-tim(1));
      fprintf(' seconds \n');
    end
  end
end

quit; % needed to run automatic updates with full desktop mode
      % in order to get image processing to work correctly.
      % added by Anton 7/21/2011

% ------------------------------------------------- Other programs to run

% break

% Turned off for automated updates - Anton 6/3/2011

% zListNonRedundantSet              % 

% zFindExemplars                     % 

% zWriteHTMLFileList
% Turned off for automated updates - Anton 6/3/2011