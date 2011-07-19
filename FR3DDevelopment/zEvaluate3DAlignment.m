% zEvaluate3DAlignment(Alignment,File) superimposes elements of Alignment
% from files File(1) and File(2)

function [void] = zEvaluate3DAlignment(Alignment,SortByMaxDiscrep,File)

if nargin < 3,
  File = zAddNTData({Alignment.Filename(1), Alignment.Filename(2)},2);
end

NumCorr = length(Alignment.Correspondence);

stop         = 0;
w            = 16;                          % row of alignment to start at
superpositionrange = 1;
neighborhood = 1;
sugar        = 1;
fontsize     = 8;
labelbases   = 8;
VP.Write     = 0;
inter        = 1;

if SortByMaxDiscrep == 0,

  CorrList = 1:NumCorr;

else

  Shift = [];

  % ------------------------------- loop through elements of alignment ------

  Discrepancy = zeros(1,NumCorr);
  VP.Plot  = 0;
  VP.Write = 0;

  for Corr = 1:NumCorr,                            

    % ------------ Make lists of nucleotides shared by the molecules

    startCorr = max(1,Corr-neighborhood);
    endCorr   = min(NumCorr, Corr + neighborhood);

    Indices = [];

    for r = startCorr:endCorr,
      Indices = [Indices [Alignment.Correspondence(r).File1; Alignment.Correspondence(r).File2]];
    end

    Distinct = ones(1,length(Indices(1,:)));
    for k = 1:length(Distinct),
      for m = 1:(k-1),
         if all(Indices(:,m) == Indices(:,k)),
           Distinct(k) = 0;
         end
      end
    end

    Indices = Indices(:,find(Distinct));

    L = length(Indices(1,:));                  % number being superimposed

    % ------------ Superimpose nucleotides
    if length(Indices(1,:)) > 2,
      [d,s] = zSuperimposeNucleotides(File(1),Indices(1,:),File(2),Indices(2,:),VP,L);
      Discrepancy(Corr) = d;
%      Shift = [Shift; s];
    end

    % ------------ Display nucleotide lists being superimposed

    for j = 1:2,
%      fprintf('%s nucleotides', Alignment.Filename{j});
      for m = 1:length(Indices(j,:)),
%        fprintf(' %s', File(j).NT(Indices(j,m)).Number);
      end
%      fprintf('\n');
    end

    if length(Indices(1,:)) > 2,
%      fprintf('Geometric discrepancy %8.4f\n\n', d);
    else
%      fprintf('Not enough nucleotides to superimpose\n\n');
    end
  end

  [d,CorrList] = sort(-1*SortByMaxDiscrep*Discrepancy);         % put highest discrepancy first

hist(Discrepancy)

size(CorrList)

end

% ------------------------------- make note of correspondences indicated

AMatrix = sparse(zeros(File(1).NumNT,File(2).NumNT));

for i = 1:NumCorr,
  for j = 1:length(Alignment.Correspondence(i).File1),
    AMatrix(Alignment.Correspondence(i).File1(j),Alignment.Correspondence(i).File2(j)) = 1;
  end
end


% ------------------------------- display menu -----------------------------

VP.Plot = 1;                          % plot the bases
stop = 0;
w = 1;                                % starting element of the alignment

while stop == 0,                            

  k=menu('3D alignment','Next element','Previous element', ... % 1,2
         'Change Superposition Range', ...                     % 3
         'Change Additional Neighborhood', ...                 % 4
         'Toggle sugar','Toggle display', ...                  % 5,6
         'Write to PDB', ...                                   % 7
         'Show Interaction Matrices', 'Quit');                 % 8,9

  switch k                               % k is the menu choice
    case 1                                      % next row of alignment
      w = w + 1;

    case 2                                      % previous row of alignment
      w = w - 1;
      if w < 1, 
        w = NumCorr; 
      end

    case 3                                      % enlarge neighborhood
      sr = [2 3 4 1];                           % neighborhood setting list
      superpositionrange = sr(superpositionrange);

    case 4
      neigh = [1 2 3 0];                        % neighborhood setting list
      neighborhood  = neigh(1 + neighborhood);

    case 5                                      % toggle sugar
      sugar = 1 - sugar;

    case 6                                      % toggle numbers
      if labelbases == 0,
        labelbases = fontsize;
      elseif labelbases > 0,
        labelbases = 0;
      end

    case 7                                      % write to PDB
      VP.Write = 1;

    case 8                                      % show interaction matrices
      inter = 1 - inter;

    case 9                                      % quit
      stop = 1;

  end  % switch statement for menu

  VP.Sugar = sugar;
  VP.LabelBases = labelbases;

  VP.LineStyle     = '-';
  VP.LineThickness = '1';
  VP.ConnectSugar  = 1;
  VP.Grid          = 0;

  if any([1 2 3 4 5 6] == k),

    % ------------ Make lists of nucleotides shared by all molecules

    startCorr = max(1,CorrList(w)-superpositionrange);
    endCorr   = min(NumCorr, CorrList(w) + superpositionrange);

    Indices = [];

    for r = startCorr:endCorr,
      Indices = [Indices [Alignment.Correspondence(r).File1; Alignment.Correspondence(r).File2]];
    end

    Distinct = ones(1,length(Indices(1,:)));

    for k = 1:length(Distinct),
      for m = 1:(k-1),
         if all(Indices(:,m) == Indices(:,k)),
           Distinct(k) = 0;
         end
      end
    end

    Indices = Indices(:,find(Distinct));

    L = length(Indices(1,:));                  % number being superimposed

    E1 = zExpandList(Indices(1,:),neighborhood,File(1).NumNT);
    E2 = zExpandList(Indices(2,:),neighborhood,File(2).NumNT);

    for e = 1:length(E1),
      if sum(AMatrix(E1(e),E2)) == 0,
        fprintf('Nucleotide %s%s in %s has no correspondence listed\n', File(1).NT(E1(e)).Base, File(1).NT(E1(e)).Number, File(1).Filename);
      end
    end

    for e = 1:length(E2),
      if sum(AMatrix(E1,E2(e))) == 0,
        fprintf('Nucleotide %s%s in %s has no correspondence listed\n', File(2).NT(E2(e)).Base, File(2).NT(E2(e)).Number, File(2).Filename);
      end
    end


    Indices1 = [Indices(1,:) E1];
    Indices2 = [Indices(2,:) E2];

    if inter == 0,

      % ------------ Display nucleotide lists being superimposed

      for j = 1:2,
        fprintf('%s nucleotides', Alignment.Filename{j});
        for m = 1:length(Indices(j,:)),
          fprintf(' %s', File(j).NT(Indices(j,m)).Number);
        end
        fprintf('\n');
      end

    else

      % ------------ Show interaction matrices

      for k=1:2,
        zShowInteractionTable(File(k),Indices(k,:));
      end

    end

    % ------------ Superimpose and display nucleotides

    if length(Indices(1,:)) > 2,
      figure(1);
      [az,el] = view;
      clf
      [d,s] = zSuperimposeNucleotides(File(1),Indices1,File(2),Indices2,VP,L);
      view(az,el);
      fprintf('Geometric discrepancy %8.4f\n\n', d);
    else
      fprintf('Increase size of neighborhood\n');
    end

  end

  rotate3d on
  drawnow

  VP.Write = 0;

end  % end while



