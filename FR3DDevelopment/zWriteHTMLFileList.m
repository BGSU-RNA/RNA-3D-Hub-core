% zWriteHTMLFileList produces large HTML tables of all PDB files with RNA and with only the non-redundant dataset

function [void] = zWriteHTMLFileList(reportdate)

RepAndBetter = 0;               % if > 1, only list representative and better structures

webroot = 'http://rna.bgsu.edu/FR3D/AnalyzedStructures'; % Anton 15/2/2011

if nargin < 1,
  reportdate = datestr(date,'yyyy-mm-dd');
else
  reportdate = reportdate;
end

destination = ['Web' filesep 'AnalyzedStructures']; % Anton 1/18/2011 To make mac-compatible
if ~exist(destination,'dir')  % Anton 1/18/2011 To make mac-compatible
	mkdir(destination);
end

load PDBInfo
% load([pwd filesep 'FR3DSource' filesep 'PDBInfo.mat']);

for r = 1:length(t(:,1)),
  if isempty(t{r,10}),
    t{r,10} = t{r,1};        % in case there is no representative for this file
  end
end

Reports = [[1 Inf]; [2 Inf]; [2 20]; [2 4]; [2 3.5]; [2 3]; [2 2.5]; [2 2]; [2 1.5]; [2 1]];

for r = 1:length(Reports(:,1)),

  fprintf('Writing the html file list up to resolution %7.1f\n', Reports(r,2));

  Format = Reports(r,1);
  MaxRes = Reports(r,2);

  % Format = 1 writes the page of all RNA-containing PDB files
  % Format = 2 writes the page of non-redundant PDB files

  % =========================== Write headers of web pages

  switch Format
  case 1,                                             % list all files
    i = find(n(:,2) > 0);                             % files with nucleotides
    i = 1:length(t(:,1));                             % list all structures
    [tt,nn] = zEquivalents(t(i,:),n(i,:),5,0);        % files equivalent to these
    fid = fopen([destination filesep 'index.html'],'w');    
    fprintf(fid,'<html>');
    fprintf(fid,'<head>\n');
    fprintf(fid,'<title>List of all %d PDB structures containing at least one RNA nucleotide as of %s</title>\n',length(tt(:,1)),reportdate);
    fprintf(fid,'<link rel="stylesheet" media="all" type="text/css" href="http://rna.bgsu.edu/nrlist/css/nrlists.css">\n');
    fprintf(fid,'<script src="http://rna.bgsu.edu/nrlist/js/sorttable.js"></script>\n');  %Anton 5/20/11   
    fprintf(fid,'</head>\n');
    fprintf(fid,'<body>\n');
    fprintf(fid,'<h3>List of all %d RNA-containing structures as of %s</h3>\n',length(tt(:,1)),reportdate);
    fprintf(fid,'Equivalent structures are grouped together, with the representative listed first and the others in decreasing order of basepairs per nucleotide.<br>\n');
    fprintf(fid,'<a href="All_structures_%s_list.pdb">Click here for the list file to use with FR3D</a> (Save in FR3D/PDBFiles.)<br>', reportdate);
    fprintf(fid,'Click here for the non-redundant list of high-resolution x-ray structures up to: ');
    for a = length(Reports(:,1)):-1:2,
      MR = [num2str(Reports(a,2)) 'A'];
      MRText = strrep(num2str(MR),'.',',');
      if Reports(a,2) == Inf,
        MR     = 'All Resolutions';
        MRText = 'All_Resolution';
      end
      fprintf(fid,'%s',['<a href="Nonredundant_' MRText '.html">' MR '</a> ']);
    end
    fprintf(fid,'\n');

    % ------------ Study: Only write the representative and structures with better resolution

    if RepAndBetter > 1,

      i = find(nn(:,1)<=4);                          % 4A or better
      tt = tt(i,:);
      nn = nn(i,:);

      currentres = nn(1,1);                          % res of first structure
      currentrep = 1;                                % index of rep
      Keep = [];
      for j = 2:length(tt(:,1)),                     % go through structures
        
        if strcmp(tt{j,10},tt{j-1,10}),              % same class as previous
          if nn(j,1) < floor(2*currentres-0.001)/2,  % much better resolution
%          if nn(j,1) < currentres,                   % better resolution
            Keep(j) = 1;                             % keep this one
            Keep(currentrep) = 1;                    % keep rep too
          end
        end

        if ~strcmp(tt{j,10},tt{j-1,10}),
          currentres = nn(j,1);
          currentrep = j;
        end
      end
      i = find(Keep);
      tt = tt(i,:);
      nn = nn(i,:);
    end



  case 2,
    i = find((n(:,1) > 0) .* (n(:,1) <= MaxRes) .* (n(:,2) > 0));   
                           % Not NMR, res to MaxRes, at least one nucleotide

    if MaxRes == Inf,           
      i = find((n(:,2) > 0));   % allow NMR, have at least one basepair
    end

    [tt,nn] = zEquivalents(t(i,:),n(i,:),5,1);      % get equivalents

    Keep = [];
    Keep(1) = 1;
    for j = 2:length(tt(:,1)),
      if ~strcmp(tt{j,10},tt{j-1,10}),              % first instance from group
        Keep(j) = 1;
      end
    end

    i = find(Keep);
    tt = tt(i,:);
    nn = nn(i,:);

    MR = [num2str(Reports(r,2)) 'A'];
    MRText = strrep(num2str(MR),'.',',');
    if MaxRes == Inf,
      MRText = 'All_Resolution';
      MR     = 'All Resolutions';
    end

    fid = fopen([destination filesep 'Nonredundant_' MRText '.html'],'w');
    fprintf(fid,'<html>');
    fprintf(fid,'<head>\n');
    fprintf(fid,'<link rel="stylesheet" media="all" type="text/css" href="http://rna.bgsu.edu/nrlist/css/nrlists.css">\n');
    fprintf(fid,'<script src="http://rna.bgsu.edu/nrlist/js/sorttable.js"></script>\n');  %Anton 5/20/11        
    if MaxRes < Inf,
      fprintf(fid,'<title>Non-redundant list of PDB structures from X-ray crystallography with resolution %7.1fA or better as of %s</title>\n', MaxRes, reportdate);
    else
      fprintf(fid,'<title>Non-redundant list of PDB structures of all resolutions as of %s</title>\n',reportdate);
    end
    fprintf(fid,'</head>\n');
    fprintf(fid,'<body>\n');
    if MaxRes < Inf,
      fprintf(fid,'<h3>Non-redundant list of %d PDB structures from X-ray crystallography with resolution &le; %7.1fA and containing at least one RNA nucleotide, as of %s</h3>\n', length(i), MaxRes, reportdate);
    else
      fprintf(fid,'<h3>Non-redundant list of %d PDB structures of all resolutions and containing at least one RNA nucleotide, as of %s</h3>\n', length(i),reportdate);
    end
    fprintf(fid,'They are sorted by the length of the longest unique chain, longest to shortest.<br>\n');
    fprintf(fid,'<a href="Nonredundant_%s_%s_list.pdb">Click here for the list file to use with FR3D</a> (Save in FR3D/PDBFiles.)<br>', MRText, reportdate);
    fprintf(fid,'<a href="index.html">Click here for the complete list</a><br>\n');
    fprintf(fid,'Click here for the non-redundant list of high-resolution x-ray structures up to: ');
    for a = length(Reports(:,1)):-1:2,
      MR = [num2str(Reports(a,2)) 'A'];
      MRTextt = strrep(num2str(MR),'.',',');
      if Reports(a,2) == Inf,
        MR     = 'All Resolutions';
        MRTextt = 'All_Resolution';
      end
      fprintf(fid,'%s',['<a href="Nonredundant_' MRTextt '.html">' MR '</a> ']);
    end
    fprintf(fid,'\n');

  end

  % ================================= Write the table header

  fprintf(fid,'<table border=1 class="sortable">'); %Anton 5/20/11

  fprintf(fid,'<tr><th class="narrow">PDBID</th><th class="narrow">FR3D analysis</th><th class="narrow">Chain(s)</th><th class="narrow">NTs</th><th class="narrow">Pairs</th><th class="narrow">Res</th><th class="narrow">Date</th><th>Organism</th><th>Description provided by author</th><th>Author</th>');

  switch Format,
  case 1,
    fprintf(fid,'<th class="narrow">Represented by</th></tr>\n');
  case 2,
    fprintf(fid,'<th>Also represents</th></tr>\n');
  end

  % ================================= Write the table entries

  for r = 1:length(tt(:,1)),
    fprintf(fid,'<tr>');
    fprintf(fid,'<td><a href="http://www.rcsb.org/pdb/cgi/explore.cgi?pdbId=%s" title="Link to PDB">%s</a></td>', upper(tt{r,1}), upper(tt{r,1}));

    fprintf(fid,'<td><a href="%s/%s/index.html" title="Link to FR3D analysis. %0.4f pairs per nucleotide">%s</a></td>', webroot,upper(tt{r,1}),nn(r,3)/nn(r,2),tt{r,1}); % CLZ 2011-05-17

    if nn(r,1) == 0 || isnan(nn(r,1)),
      res = tt{r,3};
    else
      res = sprintf('%0.1f',nn(r,1));
    end
%  res = strrep(res,'ELECTRON MICROSCOPY','EM');
%  res = strrep(res,'SOLUTION NMR','NMR');

    fprintf(fid,'<td align="center">%s</td>', char(tt{r,12}));
    fprintf(fid,'<td align="center">%d</td><td align="center">%d</td><td align="center">%s</td>', nn(r,2), nn(r,3), res);

    org = tt{r,8};
    if isempty(org),
      org = '&nbsp;';
    end

    [year month day hh mm ss] = datevec(tt{r,4}, 'dd-mm-yy'); % Anton 2011-05-24, to make dates sortable
    sortabledate = sprintf('%i%02d%02d%02d%02d%02d',year,month,day,hh,mm,ss); % Anton 2011-05-24, to make dates sortable

    fprintf(fid,'<td sorttable_customkey="%s">%s</td><td>%s</td><td>%s</td><td>%s</td>',sortabledate,tt{r,4},org,tt{r,2},tt{r,5}); % Anton 2011-05-24, to make dates sortable
%     fprintf(fid,'<td>%s</td><td>%s</td><td>%s</td><td>%s</td>',tt{r,4},org,tt{r,2},tt{r,5}); % Anton 2011-05-24, to make dates sortable

    switch Format,
    case 1,
      fprintf(fid,'<td>');
      if ~isempty(tt{r,10}),
        fprintf(fid,'<a href="%s/%s/index.html">%s</a> ', webroot,upper(tt{r,10}),tt{r,10}); % CLZ 2011-05-17
      else
        fprintf(fid,'&nbsp;');
      end
      fprintf(fid,'</td>');
    case 2,
      g = find(ismember(t(:,10),tt{r,10}));         % find *all* others in this group, regardless of resolution
      fprintf(fid,'<td>');
      if length(g) > 1,
        [y,i] = sort(-n(g,3)./n(g,2));                       % sort by pairs per nucleotide
        g = g(i);
        for h = 1:length(g),
          if ~strcmpi(t{g(h),1},tt{r,1}),
            if n(g(h),1) == 0 || isnan(n(g(h),1)),
              res = t{g(h),3};
            else
              res = sprintf('Resolution %0.1fA',n(g(h),1));
            end
            fprintf(fid,'<a href="%s/%s/index.html" title="%s, %d nucleotides, %d basepairs, %0.4f pairs per nucleotide">%s</a> ', webroot,upper(t{g(h),1}),res,n(g(h),2),n(g(h),3),n(g(h),3)/n(g(h),2),t{g(h),1}); % CLZ 2011-05-17
          end
        end

        fprintf(fid,'<a href="%s/%s/%s_compare_circular.html">Compare</a>', webroot, tt{r,1}, tt{r,1});

      else
        fprintf(fid,'&nbsp;');
      end
      fprintf(fid,'</td>');
    end

    fprintf(fid,'</tr>\n');
  end

  fprintf(fid,'</table>\n');
  fprintf(fid,'</html>');

  fclose(fid);

  % ====================================== Write NR lists for use in FR3D

  switch Format
  case 1,
    fid2 = fopen([destination filesep 'All_structures_' reportdate '_list.pdb'],'w');
    for i = 1:length(tt(:,1)),
      fprintf(fid2,'%s\n', tt{i,1});                % write to the list
    end
    fclose(fid2);

  case 2,
    if MaxRes == Inf,
      MRText = 'All_Resolution';
    end

    fid2 = fopen([destination filesep 'Nonredundant_' MRText '_' reportdate '_list.pdb'],'w');

    for i = 1:length(tt(:,1)),
      fprintf(fid2,'%s\n', tt{i,1});                % write to the list
    end
    fclose(fid2);

  end
end
