fprintf('Superimposing structures with sequence similarity.\n',100*p);

done = zeros(1,F);                % whether each file has been considered
Discrepancies = sparse(F,F);

stim = cputime;

for i = 1:(F-1),                  % loop through files
  if done(i) == 0,                % if file i hasn't been considered yet
    j = find(Close(i,:));         % files that are "close" to i
    if length(j) < 2,             % no file is close to i
      done(i) = 1;                % file i has a unique sequence
    else

%      File = zAddNTData(t(j,1));     % load nucleotide data
                                     % extended by transitivity so that,
                                     % hopefully, each file only needs to be
                                     % loaded once

      clear File

      for jj = 1:length(j),
        FF = zAddNTData(t{j(jj),1},0,[],1);
        FFF.NT = FF.NT;
        FFF.NumNT = FF.NumNT;
        FFF.Filename = FF.Filename;
        FFF.Info = FF.Info;
        FFF.LongestChain = FF.LongestChain;
        File(jj) = FFF;
      end

      for m = 1:length(j),
        fprintf('%4s has %4d nucleotides, %4d in longest chain, %4d basepairs, ', File(m).Filename, File(m).NumNT, n(j(m),4), n(j(m),3));
        if isempty(File(m).Info.Resolution),
          fprintf('resolution  ---- ');
        else
          fprintf('resolution %5.2f ', File(m).Info.Resolution);
        end

        Info = File(m).Info;

        fprintf(' %10s | %s | %s | %s\n', Info.ReleaseDate, Info.Source, Info.Descriptor, Info.Author);
      end

      fprintf('Percent agreement of base sequence, using alignment of sequences\n');

      fprintf('           ');
      for m = 1:length(j),
        fprintf(' %4s  ', File(m).Filename);
      end
      fprintf('\n');

      for m = 1:length(j),
        done(j(m)) = 1;
        T = sprintf('%4s', File(m).Filename);
        if isempty(File(m).Info.Resolution),
          T = [T sprintf(' ---- ')];
        else
          T = [T sprintf('%5.2f ', File(m).Info.Resolution)];
        end
        for nn = 1:length(j),
          if nn == m || prop(j(m),j(nn)) == 0,
            T = [T sprintf('       ')];
          else
            T = [T sprintf(' %5.1f ', 100*full(prop(j(m),j(nn))))];
          end
        end
        fprintf('%s\n',T);
      end
      fprintf('\n');

      fprintf('Sequences of longest chains:\n');
      for m = 1:length(j),                        % print sequences
        if length(t{j(m),11}) < 10000,
          fprintf('%4s %s\n', File(m).Filename, t{j(m),11});
        end
      end
      fprintf('\n');

      discmat = zeros(length(j));

if length(j) > 30,
  fprintf('Calculating discrepancies\n');
  drawnow
end

      for m = 1:length(j),
        for nn = (m+1) : length(j),
          if ~isempty(align{j(m),j(nn)}),
            malign = align{j(m),j(nn)}(1,:);     % use stored data
            nalign = align{j(m),j(nn)}(2,:);
          else
            malign = [];
            nalign = [];
          end

          if isempty(malign) || closeseq(j(m),j(nn)) == 0,
            d = Inf;
          else
            d = xDiscrepancy(File(m),File(m).LongestChain(1)-1+malign,File(nn),File(nn).LongestChain(1)-1+nalign);
          end

          discmat(m,nn) = d;
          Discrepancies(j(m),j(nn)) = d;
          Discrepancies(j(nn),j(m)) = d;

          if ~(d <= maxd),                  % allow for d = NaN
            closeseq(j(m),j(nn)) = 0;       % these are not that close!
            closeseq(j(nn),j(m)) = 0;  
          end     
        end
      end

      discmat = discmat + discmat';

      fprintf('Geometric discrepancy between aligned bases, where there is sequence similarity\n');
      fprintf('           ');

      for m = 1:length(j),
        fprintf(' %4s  ', File(m).Filename);
      end
      fprintf('\n');

      for m = 1:length(j),
        T = sprintf('%4s', File(m).Filename);
        if isempty(File(m).Info.Resolution),
          T = [T sprintf(' ---- ')];
        else
          T = [T sprintf('%5.2f ', File(m).Info.Resolution)];
        end
        for nn = 1:length(j),
          if nn == m,
            T = [T sprintf('       ')];
          elseif discmat(m,nn) < Inf,
            T = [T sprintf(' %5.2f ', discmat(m,nn))];
          else
            T = [T sprintf('       ')];
          end
        end
        fprintf('%s\n',T);
      end
      fprintf('\n');
    end      
  end
  drawnow
end

fprintf('\nSuperimposing structures took %7.1f minutes.\n\n', (cputime-stim)/60);
