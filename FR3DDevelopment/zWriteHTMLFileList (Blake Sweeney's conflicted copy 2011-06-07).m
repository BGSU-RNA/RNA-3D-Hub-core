% zWriteHTMLFileList produces large HTML tables of all PDB files with RNA and with only the non-redundant dataset

function [void] = zWriteHTMLFileList(reportdate)


webroot = 'http://rna.bgsu.edu/FR3D/AnalyzedStructures'; % Anton 15/2/2011

if nargin < 1,
  reportdate = '2011-01-07';
else
  reportdate = reportdate;
end

destination = ['Web' filesep 'AnalyzedStructures']; % Anton 1/18/2011 To make mac-compatible
if ~exist(destination,'dir')  % Anton 1/18/2011 To make mac-compatible
	mkdir(destination);
end

% load PDBInfo
load([pwd filesep 'FR3DSource' filesep 'PDBInfo.mat']);

for r = 1:length(t(:,1)),
  if isempty(t{r,10}),
    t{r,10} = t{r,1};
  end
end

Reports = [[2 4]; [1 Inf]; [2 Inf]; [2 20]; [2 3.5]; [2 3]; [2 2.5]; [2 2]; [2 1.5]];

for r = 1:length(Reports(:,1)),

  fprintf('Writing the html file list up to resolution %7.1f\n', Reports(r,2));

  Format = Reports(r,1);
  MaxRes = Reports(r,2);

  % Format = 1 writes the page of all RNA-containing PDB files
  % Format = 2 writes the page of non-redundant PDB files

  switch Format
  case 1, 
    i = find(n(:,2) > 0);                             % files with nucleotides
    [tt,nn] = zEquivalents(t(i,:),n(i,:),5,0);        % files equivalent to this
%    fid = fopen('Web\AnalyzedStructures\index.html','w');
    fid = fopen([destination filesep 'index.html'],'w');    
    fprintf(fid,'<html>');
    fprintf(fid,'<head>\n');
    fprintf(fid,'<title>List of all PDB structures containing RNA nucleotides as of %s</title>\n',reportdate);
    fprintf(fid,'</head>\n');
    fprintf(fid,'<body>\n');
    fprintf(fid,'<h3>List of all RNA-containing structures as of %s</h3>\n',reportdate);
    fprintf(fid,'Equivalent structures are grouped together.<br>\n');
    fprintf(fid,'<a href="All_structures_%s_list.pdb">Click here for the list file to use with FR3D</a> (Save in FR3D/PDBFiles.)<br>', reportdate);
    fprintf(fid,'Click here for the non-redundant list of high-resolution x-ray structures up to: ');
    for a = 2:length(Reports(:,1)),
      MR = [num2str(Reports(a,2)) 'A'];
      MRText = strrep(num2str(MR),'.',',');
      if Reports(a,2) == Inf,
        MR     = 'All Resolutions';
        MRText = 'All_Resolution';
      end
      fprintf(fid,'%s',['<a href="Nonredundant_' MRText '.html">' MR '</a> ']);
    end
    fprintf(fid,'\n');

  case 2,
    i = find((n(:,1) > 0) .* (n(:,1) <= MaxRes) .* (n(:,3) > 0));   
                                % Not NMR, res to MaxRes, at least one pair

    if MaxRes == Inf,           
      i = find((n(:,3) > 0));   % allow NMR, have at least one basepair
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
    if MaxRes < Inf,
      fprintf(fid,'<title>Non-redundant list of PDB structures from X-ray crystallography with resolution %7.1fA or better as of %s</title>\n', MaxRes, reportdate);
    else
      fprintf(fid,'<title>Non-redundant list of PDB structures of all resolutions as of %s</title>\n',reportdate);
    end
    fprintf(fid,'</head>\n');
    fprintf(fid,'<body>\n');
    if MaxRes < Inf,
      fprintf(fid,'<h3>Non-redundant list of %d PDB structures from X-ray crystallography with resolution &le; %7.1fA and containing at least one RNA basepair, as of %s</h3>\n', length(i), MaxRes, reportdate);
    else
      fprintf(fid,'<h3>Non-redundant list of %d PDB structures of all resolutions and containing at least one RNA basepair, as of %s</h3>\n', length(i),reportdate);
    end
    fprintf(fid,'They are sorted by the length of the longest unique chain, longest to shortest.<br>\n');
    fprintf(fid,'<a href="Nonredundant_%s_%s_list.pdb">Click here for the list file to use with FR3D</a> (Save in FR3D/PDBFiles.)<br>', MRText, reportdate);
    fprintf(fid,'<a href="index.html">Click here for the complete list</a>\n');
  end

  fprintf(fid,'<table border=1>');

  fprintf(fid,'<tr><th>PDBID</th><th>FR3D analysis</th><th>NTs</th><th>Pairs</th><th>Res</th><th>Date</th><th>Organism</th><th>Description provided by author</th><th>Author</th>');

  switch Format,
  case 1,
    fprintf(fid,'<th>Represented by</th></tr>\n');
  case 2,
    fprintf(fid,'<th>Also represents</th></tr>\n');
  end

  for r = 1:length(tt(:,1)),
    fprintf(fid,'<tr>');
    fprintf(fid,'<td><a href="http://www.rcsb.org/pdb/cgi/explore.cgi?pdbId=%s">%s</a></td>', upper(tt{r,1}), upper(tt{r,1}));
    fprintf(fid,'<td><a href="%s/%s/index.html">FR3D</a></td>', webroot,upper(tt{r,1})); % Anton 15/2/2011

    if nn(r,1) == 0 || isnan(nn(r,1)),
      res = tt{r,3};
    else
      res = sprintf('%0.1f',nn(r,1));
    end
%  res = strrep(res,'ELECTRON MICROSCOPY','EM');
%  res = strrep(res,'SOLUTION NMR','NMR');

    fprintf(fid,'<td align="center">%d</td><td align="center">%d</td><td align="center">%s</td>', nn(r,2), nn(r,3), res);

    org = tt{r,8};
    if isempty(org),
      org = '&nbsp;';
    end

    fprintf(fid,'<td>%s</td><td>%s</td><td>%s</td><td>%s</td>',tt{r,4},org,tt{r,2},tt{r,5});

    switch Format,
    case 1,
      fprintf(fid,'<td>');
      if strcmpi(tt{r,10},tt{r,1}),
        fprintf(fid,'%s',tt{r,1});
      elseif ~isempty(tt{r,10}),
       fprintf(fid,'%s',tt{r,10});
      else
        fprintf(fid,'&nbsp;');
      end
      fprintf(fid,'</td>');
    case 2,
      g = find(ismember(t(:,10),tt{r,10}));         % find others in this group
%      fprintf(fid,'<td width="250">');
      fprintf(fid,'<td>');
      if length(g) > 1,
        for h = 1:length(g),
          if ~strcmpi(t{g(h),1},tt{r,1}),
            fprintf(fid,'%s ',t{g(h),1});
          end
       end
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
