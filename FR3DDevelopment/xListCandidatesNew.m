% xListCandidates(Search) prints a candidate list to the screen
% The optional argument NumToOutput limits the list's length
% The optional argument WheretoOutput has this effect:
%   Value 1 : prints a wide listing to the Matlab command window
%   Value 2 : prints a wide listing to an Editbox
%   Value 3 : prints a narrow listing to an Editbox
%   Value 5 : returns a wide listing, doesn't print anything
%   Value 6 : returns the text of the narrow listing
%   Value 7 : narrow listing with information on the organism

% The PC compiled version does both 2 and 3.

% It may be run directly from Matlab using the command:
%   xListCandidates(Search);

function [Text] = xListCandidates(Search,NumToOutput,WheretoOutput,Param)

File        = Search.File;
Candidates  = Search.Candidates;

if ~isfield(Search,'Query'),
  Query.Geometric = 0;
  Query.Name = '';
else
  Query = Search.Query;
end

[s,t]       = size(Candidates);
N           = t-1;

if s == 0,
  fprintf('There are no candidates to list\n');
  return
end

if N == 2,
  CP = zeros(1,s);
end

if nargin < 2,
  NumToOutput = Inf;                    % limit on number printed to screen
end

if (nargin < 3),
  if isdeployed,
    WheretoOutput = 2;
    xListCandidates(Search,NumToOutput,3);
  else
    WheretoOutput = 1;                    % make a wide listing
  end
end

% -------------------------------------- collect descriptor lines

if isfield(Search,'SaveName'),
  Text{1,1} = Search.SaveName;
else
  Text{1,1} = '';
end

if isfield(Search.Query,'Name'),
  Text{2,1} = Search.Query.Name;
else
  Text{2,1} = '';
end

if isfield(Search.Query,'Description'),
  Text{3,1} = Search.Query.Description;
else
  Text{3,1} = '';
end

% -------------------------------------- collect header line

t = 4;

Text{t,1} = 'Filename';

if isfield(Search,'AvgDisc'),
  Text{t,2} = 'Avg Discrep';
elseif Query.Geometric > 0,
  Text{t,2} = 'Discrepancy')];
else
  Text{t,2} = 'Number');
end

c = 2;                                     % column number

for i=1:N,                                 % one for each nucleotide
  c = c + 1;
  Text{t,c} = sprintf('%d', i);
end

c = c + 1;
Text{t,c} = 'Chains';

if isfield(Search,'GroupLabel'),
  c = c + 1;
  Text{t,c} = 'Group';
end

if any(WheretoOutput == [1 2 5]),
  for i=1:N,
    for j=(i+1):N,
      c = c + 1;
      Text{t,c} = sprintf('%s', [num2str(i) '-' num2str(j)]);
    end
  end
  
  c = c + 1;
  Text{t,c} = 'Configuration';
  
  for i=1:N,
    for j=(i+1):N,
      c = c + 1;
      Text{t,c} = sprintf('%s', [num2str(i) '-' num2str(j)]);
    end
  end

  for i=1:N,
    for j=1:N,
      c = c + 1;
      Text{t,c} = sprintf('%s', [num2str(i) '-' num2str(j)]);
    end
  end

  for i=1:N,
    for j=1:N,
      c = c + 1;
      Text{t,c} = sprintf('%s', [num2str(i) '-' num2str(j)]);
    end
  end

  for i=1:N,
    for j=(i+1):N,
      c = c + 1;
      Text{t,c} = sprintf('%s', [num2str(i) '-' num2str(j)]);
    end
  end
  
end   

if N == 2,
  Text{t} = [Text{t} sprintf(' Pair data')];
end
  
% -------------------------------------- list candidates

Config = {'A' , 'S'};

for i=1:min(s,NumToOutput),

  f = double(Candidates(i,N+1));               % file number for this candidate
  Indices = Candidates(i,1:N);                 % indices of nucleotides

  Text{i+t,1} = File(f).Filename;

  if isfield(Search,'DisttoCenter'),
    Text{i+t,2} = sprintf('%12.4f',Search.DisttoCenter(i));
  elseif Query.Geometric > 0,
    Text{i+t,2} = sprintf('%12.4f',Search.Discrepancy(i));
  else
    Text{i+t,2} = sprintf('%12d',Search.Discrepancy(i));      % original candidate number
  end

  c = 2;

  for j=1:N,                                   % print nucleotide base and num
    c = c + 1;
    Text{i+t,c} = sprintf('%1s%5s',File(f).NT(Indices(j)).Base,File(f).NT(Indices(j)).Number);
  end

  c = c + 1;

  for j=1:N,
    Text{i+t,c} = cat(2,File(f).NT(Indices).Chain);
  end

  if isfield(Search,'GroupLabel'),
    GL = Search.GroupLabel{i};
    c = c + 1;
    Text{i+t,c} = GL(1:10);
  end

  if any(WheretoOutput == [1 2 5]),
    for k=1:length(Indices),                    % list pairs and stacks
      for j=(k+1):length(Indices),
        C1 = File(f).NT(Indices(k)).Code;
        C2 = File(f).NT(Indices(j)).Code;
        c = c + 1;
        Text{i+t,c} = zEdgeText(File(f).Edge(Indices(k),Indices(j)),1,C1,C2);
      end
    end
    
    conf = '';
    for k=1:length(Indices),
      conf = [conf Config{File(f).NT(Indices(k)).Syn+1}];
    end
    c = c + 1;
    Text{i+t,c} = conf;
    
    for k=1:length(Indices),
      for j=(k+1):length(Indices),
        c = c + 1;
        Text{i+t,c} = sprintf('%6d', abs(double(Indices(k))-double(Indices(j))));
      end
    end

    for k=1:length(Indices),
      for j=1:length(Indices),
        c = c + 1;
        Text{i+t,c} = zBasePhosphateText(File(f).BasePhosphate(Indices(k),Indices(j)));
        end
      end
    end

    for k=1:length(Indices),
      for j=1:length(Indices),
        c = c + 1;
        Text{i+t,c} = zBaseRiboseText(File(f).BaseRibose(Indices(k),Indices(j)));
        end
      end
    end

    for k=1:length(Indices),
      for j=(k+1):length(Indices),
        bbc = max(File(f).Backbone(Indices(j),Indices(k)),File(f).Backbone(Indices(k),Indices(j)));
        c = c + 1;
        Text{i+t,c} = zBackboneText(bbc));
      end
    end

  end
    
  if N == 2,                        % special treatment for pairs

    CP(i) = norm(File(f).NT(Candidates(i,1)).Sugar(1,:) - ...
                          File(f).NT(Candidates(i,2)).Sugar(1,:));
    c = c + 1;
    Text{i+t,c} = sprintf('   C1*-C1*: %8.4f', CP(i));
    NT1 = File(f).NT(Candidates(i,1));
    NT2 = File(f).NT(Candidates(i,2));
    Edge= full(File(f).Edge(Candidates(i,1),Candidates(i,2)));
    c = c + 1;
    Text{i+t,c} = sprintf('%7.1f ', Edge);

    if isfield(File,'Crossing'),
      ii = Candidates(i,1);
      jj = Candidates(i,2);
      if (File(f).Crossing(ii,jj) == 0) % && abs(File(f).Edge(ii,jj)) < 15,
        r = 'Nested';
      elseif File(f).Range(ii,jj) > 0,
        r = sprintf('Crossing %4d', full(File(f).Range(ii,jj)));
      else
        r = '';
      end
      c = c + 1;
      Text{i+t,c} = r;
    end
  end

  % for pairs, add sequence variant information if available

  if isfield(File,'Nucl') && (WheretoOutput < 4),
    a = {};
    for j = 1:N,
      if ~isempty(File(f).Nucl(Candidates(i,j)).Motif),
        a = [a File(f).Nucl(Candidates(i,j)).Motif(1).Name];
      end
    end
%    u = unique(a);
    u = a;
    for uu = 1:length(u),
      c = c + 1;
      Text{i+t,c} = u{uu};
    end
  end

  if WheretoOutput == 7,
    c = c + 1;
    Text{i+t,c} = [Text{i+t,c} ' ' File(f).Info.Source ' | ' File(f).Info.Descriptor];
  end
end

% -------------------------------------- Remove empty columns

Empty = ones(1,length(Text(1,:)));       % presume a column is empty
Empty(1) = 0;
Empty(2) = 0;                            % these can't be empty

for c = 3:length(Text(1,:)),             % go through remaining columns
  for r = 5:length(Text(:,1)),           % and all rows
    if ~isempty(Text{r,c}),
      Empty(c) = 0;
    end
  end
end

OKText = Text(:,find(Empty == 0));

% -------------------------------------- Additional notifications and info
if (Query.Geometric > 0),
  if (Query.RelCutoff > Query.DiscCutoff) && ~isfield(Search,'AvgDisc'),
    L = length(Text);
    Extra{1} = sprintf('Some motifs with discrepancy between %7.4f and %7.4f might not appear above\n\n', Query.DiscCutoff, Query.RelCutoff);
  end
end

if s > NumToOutput,
  L = length(Text);
  Extra{2} = sprintf('Only the first %d candidates were listed.\n', NumToOutput);
end

if (N == 2) && (WheretoOutput < 4) && 0 > 1, 
  figure(12)
  clf
  hist(CP,30)
  title('Histogram of C1''-C1'' distance for these pairs');
  fprintf('Average C1''-C1'' distance is: %8.4f\n', mean(CP));
end





% -------------------------------------- Display the listing

if WheretoOutput == 3,
  mEditbox(Text,'List of Candidates',10);
elseif WheretoOutput == 2,
  mEditbox(Text,'Wide list of Candidates',7);
elseif any(WheretoOutput == [1 7]),
  for i=1:length(Text),
    fprintf('%s\n',Text{i});
  end
end
