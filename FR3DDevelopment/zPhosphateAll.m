% zPhosphateAll checks all nearby pairs of bases for base-phosphate
% interactions and stores them in a sparse matrix field BasePhosphate


%function [File,D] = zPhosphateInteractions(File,Verbose)

%if nargin == 1,
%  Verbose = 0;
%end

t = cputime;

% ------------------------------ specify cutoffs for classification

CarbonDist    = 4.0;                           % max massive - oxygen distance
nCarbonDist   = 4.5;                           % near category

NitrogenDist  = 3.5;                           % max massive - oxygen distance
nNitrogenDist = 4.0;                           % near category

AL = 130;                                      % angle limit for BP
nAL = 110;                                     % angle limit for nBP

DL([4 6 11]) = NitrogenDist;
DL([7 8 9])  = CarbonDist;

nDL([4 6 11]) = nNitrogenDist;
nDL([7 8 9])  = nCarbonDist;

% ------------------------------- define basic data

cc = 0;                           % counts how many interactions added so far
D = zeros(100000,20);            % where to save data if Verbose
T = [];                          % data on which hydrogen with which oxygen(s)

zStandardBases
Sugar = {'C1*','C2*','O2*','C3*','O3*','C4*','O4*','C5*','O5*','P','O1P','O2P','O3 of next'};

p   = [9 13 11 12];                           % rows of the phosphate oxygens
pn  = {'O5*','O3*','O1P','O2P'};              % names of phosphate oxygens

% ------------------------------ loop through files and classify

for f = 1:length(File),

  if isempty(File(f).Distance),
    c = cat(1,File(f).NT(1:File(f).NumNT).Center); % nucleotide centers
    File(f).Distance = zMutualDistance(c,16); % compute distances < 16 Angs
  end

  File(f).BasePhosphate = sparse(zeros(File(f).NumNT));

  % -------- First screening of base pairs ----------------------------------- 

  DistCutoff = 16;                              % max distance for interaction
  [i,j] = find((File(f).Distance < DistCutoff).*(File(f).Distance > 0)); 
                                                % screen by C-C distance

  i = [i; (1:length(File(f).NT))'];             % allow self interactions
  j = [j; (1:length(File(f).NT))'];             % allow self interactions

  % -------- Screen and analyze pairs ----------------------------------------

  for k = 1:length(i),                          % loop through possible pairs
   N1 = File(f).NT(i(k));                       % nucleotide i information
   N2 = File(f).NT(j(k));                       % nucleotide j information

   ph = (N2.Sugar(10,:)-N1.Fit(1,:)) * N1.Rot;  % phosphorus displacement
   if abs(ph(3)) < 4.5,                         % phosphorus close to plane

    dis = zDistance(N1.Fit, N2.Sugar(p,:)); % distances between mass & O's

    [v,pp] = min(dis,[],2);           % indices of nearest oxygen

    Dist  = [];

    for n = 1:length(v),             % massive atom to consider
      if v(n) < 4.5,                 % oxygen is close enough
          % store information for later display
          ox = (N2.Sugar(p(pp(n)),:)-N1.Fit(1,:)) * N1.Rot; % oxygen displ

          a = [f i(k) j(k) N1.Code 0 n pp(n) 0 0 ox ph File(f).Distance(i(k),j(k)) 0 0 0 0];

          % Columns:
          % 1  file number
          % 2  index of base
          % 3  index of nucleotide using phosphate
          % 4  code of base
          % 5  classification number for this massive-oxygen pair
          % 6  which massive atom is interacting
          % 7  which oxygen is interacting
          % 8  angle of interaction, in degrees
          % 9  distance from massive atom to oxygen, in Angstroms
          %10  displacement of oxygen atom relative to C1' of base
          %13  displacement of phophorus atom relative to C1' of base
          %16  distance between centers of the two bases
          %17  1 if this is the best oxygen for this hydrogen, 0 otherwise
          %18  approximate quality of the hydrogen bond
          %19  angle made by massive, hydrogen, phosphorus
          %20  distance from massive to phosphorus

          cc = cc + 1;
          if cc > length(D(:,1))
            D = [D; zeros(10000,20)];
            fprintf('Found %d interactions so far\n', cc);
          end
          D(cc,:) = a;

      end % if oxygen is near enough to an atom
    end   % loop over massive atoms
   end    % if vertical displacement of phosphate is less than 6 Angstroms
  end     % loop over nucleotide pairs
end       % loop over files


D = D(1:cc,:);                      % remove extra rows

if Verbose > 1,

  zPhosDisplay                     % temporary arrangement

end

return

% File = zAddNTData({'1s72','1j5e','2avy','2aw4','2j01'});
% zPhosphateInteractions(File,3);

Res = [];
for f = 1:length(File);
  if isempty(File(f).Info.Resolution),
    Res(f) = 10;
  else
    Res(f) = File(f).Info.Resolution;
  end
end
File = File(find(Res <= 3.0));




i = find(T(:,4) == 12);
Search = [T(i,2) T(i,3) T(i,1)];
xDisplayCandidates(File,Search)
