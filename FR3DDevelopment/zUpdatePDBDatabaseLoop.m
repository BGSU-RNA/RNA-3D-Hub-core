
for i = current:length(t(:,1)),

  fprintf('Updating structure %d, which is %s\n', i, t{i,1});

  current = i;

  try
    File = zAddNTData(t{i,1},0,[],1);          % load file
  catch
    delete(['PrecomputedData' filesep t{i,1} '.mat']);
    File = zAddNTData(t{i,1},0,[],1);          % load file
  end

  if length(File.NT) > 1,
    c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers
    File.Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
    d = sort(nonzeros(File.Distance));

    if length(d) > 1,
      if d(min(10,length(d))) < 1,
        File = zAddNTData(t{i,1},4,[],1);     % fix NMR files
      end
    end
  end

  if ~isempty(File.NT),                      % if it has nucleotides,
    n(i,2) = length(File.NT);                % store the number of NT

    E  = abs(triu(File.Edge));
    n(i,3) = full(sum(sum((E > 0) .* (E < 13)))); % number of pairs

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

%  zSaveNTData(File);

  if mod(i,100) == 0,
    save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7
  end

end
