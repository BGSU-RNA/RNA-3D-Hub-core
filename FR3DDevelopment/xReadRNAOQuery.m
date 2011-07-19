
function [Query] = xReadRNAOQuery(Filename,Verbose)

if nargin < 2,
  Verbose = 0;
end

c = 1;

Query.Description = 'Unknown';                    % default description

fid = fopen(Filename,'r');                        % open file for reading

  L = 1;
  while L > -1
    L = fgetl(fid);
    if L > -1
      n = strfind(L,'instance_of');
      if ~isempty(n),
        t = L((n(1)+length('instance_of ')):end);    % very simple
        while t(1) == ' ',
          t = t(2:end);
        end
        s = strfind(t,' ');
        if ~isempty(s),
          Query.Description = t(1:(s-1));         % given name
        end
        L = -1;                                   % stop reading
      end
    end
  end

  L = 1;

  while L > -1,
    L = fgets(fid);
    if L > -1,
      L = lower(L);                              % convert to lowercase
      if length(L) > 1,
      if ~isempty(strfind(L(1:2),'nt')),         % this line is about an nt
        a = regexp(L,'nt[0123456789]');                     % locations of nt
        b = strfind(L,'instance_of');            % specifying base

        if length(b) > 0,                        % at least one specification
          i = xExtractNum(L(a(1):end));            % number of first nt
          Query.OKCodes{i} = [0 0 0 0];               % which codes are OK
          for j = 1:length(b),                     % loop through instance_of
            z = b(j) + length('instance_of');
            while L(z) == ' ',
              z = z + 1;
            end
            switch L(z:(z+1)),
            case {'ad'}, Query.OKCodes{i}(1) = 1;
            case {'gu'}, Query.OKCodes{i}(3) = 1;
            case {'cy'}, Query.OKCodes{i}(2) = 1;
            case {'ur'}, Query.OKCodes{i}(4) = 1;
            case {'pu'}, Query.OKCodes{i}([1 3]) = [1 1];
            case {'py'}, Query.OKCodes{i}([2 4]) = [1 1];
            end
          end
        elseif length(a) > 1,                        % maybe about two nts

          i = xExtractNum(L(a(1):end));            % number of first nt
          j = xExtractNum(L(a(2):end));            % number of second nt

          Query.Edges{max(i,j),max(i,j)} = [];      % make Edge big enough
          Query.Diff{max(i,j),max(i,j)} = [];      % make Diff big enough
          if ~isempty(strfind(L,'pairs_with_')),     % particular type
            a = strfind(L,'pairs_with') + length('pairs_with_');  % location
            b = L(a:end);
            switch b(1),
            case {'c','t'},
              Query.Edges{i,j} = [Query.Edges{i,j} ' ' b(1:3)];
            case {'w','h','s'},
              Query.Edges{i,j} = [Query.Edges{i,j} ' c' b(1:2) ' t' b(1:2)];
            end

          elseif ~isempty(strfind(L,'pairs_with')),
            Query.Edges{i,j} = [Query.Edges{i,j} ' pair'];

          elseif ~isempty(strfind(L,'stack_')),
            a = strfind(L,'stack_') + length('stack_');  % location
            b = L(a:end);
            Query.Edges{i,j} = [Query.Edges{i,j} ' s' b(1) b(4)];

          elseif ~isempty(strfind(L,'stack')),
            Query.Edges{i,j} = [Query.Edges{i,j} ' stack'];
            
          elseif ~isempty(strfind(L,'cov_conn_')),
            a = strfind(L,'cov_conn_') + length('cov_conn_');  % location
            b = L(a:end);
            if (str2num(b(1)) < str2num(b(4))),              % cov_conn_3'_5'
              Query.Diff{i,j} = [Query.Diff{i,j} ' =1 <'];
            else                                             % cov_conn_5'_3'
              Query.Diff{i,j} = [Query.Diff{i,j} ' =1 >'];
            end

          elseif ~isempty(strfind(L,'5''_to')),
            Query.Diff{i,j} = [Query.Diff{i,j} ' <'];
          elseif ~isempty(strfind(L,'3''_to')),
            Query.Diff{i,j} = [Query.Diff{i,j} ' >'];

          end
        end
      end
      end
    end
  end


fclose(fid);

Query.Name = Query.Description;

N = 1;

if isfield(Query,'OKCodes'),
  N = max(N,length(Query.OKCodes));
end

if isfield(Query,'Edges'),
  [s,t] = size(Query.Edges);
  N = max(N,s);
end

if isfield(Query,'Diff'),
  [s,t] = size(Query.Diff);
  N = max(N,s);
end

if isfield(Query,'OKCodes'),
  for i = 1:N,
    if i > length(Query.OKCodes),
      Query.OKCodes{i} = [1 1 1 1];
    elseif isempty(Query.OKCodes{i}),
      Query.OKCodes{i} = [1 1 1 1];
    end
  end
end


if Verbose > 0,

Query

if isfield(Query,'Edges'),
  Query.Edges
end

if isfield(Query,'Diff'),
  Query.Diff
end

end

return

% test with xFR3DSearch

Query = xReadRNAOQuery('RNAOQuery.txt',1);
Query.SearchFiles = {'1s72'};
Query.SearchFiles = {'Nonredundant_4A_2010-05-19_list'};
Filenames = Query.SearchFiles;
Query = xConstructQuery(Query);
UsingLibrary = 1;
xFR3DSearch
xListCandidates(Search)
xDisplayCandidates(File,Search)

      if strfind(L,' '),             
      end

      if ~isempty(strfind(L,'EXPDTA')),
        Header.Expdata = L;
      end
      if strcmp(L(1:min(4,length(L))),'ATOM'),
        fprintf(out,'%s',L);
        c = c + 1;
        L = -1;                                   % stop reading header
      end

            if (i<j) && (str2num(b(1)) < str2num(b(4))),     % cov_conn_3'_5'
              Query.Edge{j,i} = [Query.Edge{i,j} ' =1 >'];
            elseif (i>j) && (str2num(b(1)) < str2num(b(4))), % cov_conn_3'_5'
              Query.Edge{i,j} = [Query.Edge{i,j} ' =1 <'];
            elseif (i<j) && (str2num(b(1)) > str2num(b(4))), % cov_conn_3'_5'
              Query.Edge{j,i} = [Query.Edge{i,j} ' =1 <'];
            elseif (i>j) && (str2num(b(1)) > str2num(b(4))), % cov_conn_3'_5'
              Query.Edge{i,j} = [Query.Edge{i,j} ' =1 >'];
            end
