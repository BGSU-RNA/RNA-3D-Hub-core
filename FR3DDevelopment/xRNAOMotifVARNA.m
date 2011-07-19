
% download the google doc spreadsheet to this file in the FR3D folder:

[n,t,r] = xlsread('RNAOMotif\RNAO motif database.xls');

numrows = length(t(:,1));

for c = 2:numrows,
  filename = t{c,1};
  fid = fopen(['RNAOMotif' filesep filename '.def'],'w');
  T = t{c,5};                           % RNAO graph description
  fprintf(fid,'%s',T);
  fclose(fid);
end

for c = 2:numrows,
  filename = ['RNAOMotif' filesep t{c,1} '.def'];

  if ~isempty(t{c,5}),                    % has a motif description
    Query = xReadRNAOQuery(filename,1);
    Query = xConstructQuery(Query);

    if strcmp(t{c,14},'RNAO_0000098'),      % hairpin, show stacks
      zMotifSignatures(Query,1,1,1);
    elseif strcmp(t{c,14},'RNAO_0000099'),  % internal loop, basepairs, rotation
      zMotifSignatures(Query,2,1,0);
    elseif strcmp(t{c,14},'RNAO_0000095'),  % paired region
      zMotifSignatures(Query,1,1,0);
    elseif strcmp(t{c,14},'RNAO_0000097'),  % unbroken stem
      zMotifSignatures(Query,1,1,0);
    elseif num2str(t{c,14}(6:end)) < 1000,  % general classes, not motifs

    else
      zMotifSignatures(Query,1,1,1);        % unknown, show once, pairs & stacks
    end

    clear T

    try
      T = xWriteVARNA(Query,[],1);
    catch
      T{1} = '<html>\n</html>\n';
      disp('************************* trouble writing VARNA description *****');
    end
    fid = fopen(['RNAOMotif' filesep t{c,1} 'VARNA.html'],'w');
    for c = 1:length(T),
      fprintf(fid,'%s\n',T{c});
    end
    fclose(fid);
  end
end
