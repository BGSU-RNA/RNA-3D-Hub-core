% zBackboneSubstitution loads the specified PDB file, then systematically
% substitutes the specified nucleotides from the specified list

function [Motif] = zBackboneSubstitution(Motif,File,AltSugar,Move)

if 0 > 1,
  MotifFile = 'Triple_cHW_tSS_GGA';
  Move      = [1];                         % nucleotide(s) to substitute
  Clash     = [1 3];                       % where the clash occurs

  MotifFile = 'Triple_tHW_tHH_AUC';
  Move      = [1 3];                       % nucleotide(s) to substitute
  Clash     = [1 3];                       % where the clash occurs

  Motif = zAddNTData(MotifFile);

  if ~exist('File'),
    File = zAddNTData('Nonredundant_2A_2011-06-25_list');
  end
end

BigDistCutoff = 4;
SmallDistCutoff = 3.5;
MinDistCutoff = 3;
MaxTry = 200000;

Verbose = 2;

VP.Sugar = 1;
VP.AtOrigin = 1;
VP.LabelBases = 0;

F = Motif;                                 % copy of the motif

% ---------------------------------------------- drop extra atom from phosphate

F.NT(1).Sugar = F.NT(1).Sugar(1:12,:);
F.NT(2).Sugar = F.NT(2).Sugar(1:12,:);
F.NT(3).Sugar = F.NT(3).Sugar(1:12,:);

% ---------------------------------------------- check for potential clashes

S1B3 = min(min(zDistance(F.NT(1).Sugar(1:12,:),F.NT(3).Fit)));
S2B3 = min(min(zDistance(F.NT(2).Sugar(1:12,:),F.NT(3).Fit)));
S3B1 = min(min(zDistance(F.NT(3).Sugar(1:12,:),F.NT(1).Fit)));
S1S3 = min(min(zDistance(F.NT(1).Sugar(1:12,:),F.NT(3).Sugar(1:12,:))));
S2S3 = min(min(zDistance(F.NT(2).Sugar(1:12,:),F.NT(3).Sugar(1:12,:))));
S2B1 = Inf;
S2S1 = Inf;

InitDist = min([S1B3 S2B3 S3B1 S1S3 S2S3 S2B1 S2S1]);

if nargin < 4,

  Move = [];

  if S1B3 < MinDistCutoff,
    Move = 1;
  end

  if S2B3 < MinDistCutoff,
    Move = [Move 2];
  end

  if S3B1 < MinDistCutoff,
    Move = [Move 3];
  end

  if S1S3 < MinDistCutoff,
    Move = [Move 1 3];
  end

  if S2S3 < MinDistCutoff,
    Move = [Move 2 3];
  end

  Move = unique(Move);

  if ~isempty(Move),
    FF = F;
    F.NT(2).Sugar = AltSugar;                   % use original sugar for 2-3

    S2B1 = min(min(zDistance(F.NT(2).Sugar(1:12,:),F.NT(3).Fit)));
    S2S1 = min(min(zDistance(F.NT(2).Sugar(1:12,:),F.NT(1).Sugar(1:12,:))));

    MMove = [];

    if S1B3 < MinDistCutoff,
      MMove = 1;
    end

    if S2B1 < MinDistCutoff,
      MMove = [MMove 2];
    end

    if S3B1 < MinDistCutoff,
      MMove = [MMove 3];
    end

    if S1S3 < MinDistCutoff,
      MMove = [MMove 1 3];
    end

    if S2S1 < MinDistCutoff,
      MMove = [MMove 2 1];
    end

    MMove = unique(MMove);
    
    if length(MMove) < length(Move),
      Move = MMove;
      S2B3 = Inf;                                % don't worry about these
      S2S3 = Inf;
      fprintf('Using sugar from NT2-NT3 pair for NT2\n');
    end
  end
end

fprintf('%s\n', F.Filename);

% =============================================== move backbones if necessary

if length(Move) > 0,                            % something needs to be moved

fprintf('Moving backbone(s) ');
for m = Move,
  fprintf('%d ',m);
end
fprintf('\n');

if Verbose > 1,
  figure(3)
  clf
  zDisplayNT(F,[2 3 1],VP);
  T = ['Original ' Motif.Filename ' model'];
  T = strrep(T,'_','-');
  title(T);
  text(0,2,0,num2str(min([S1B3 S2B3 S3B1 S1S3 S2S3])));
  view(2)
  axis off
  drawnow
  figure(4)
  clf
  drawnow
end

att = 0;                                    % number of attempts
c = 1;
Mot{1} = F;                                 % original motif
clear Mod
Mod(1,5) = 0;
Mod(1,6) = min([S1B3 S2B3 S3B1 S1S3 S2S3 S2B1 S2S1]); % current best w/ same configs

OrigConfig = cat(1,F.NT.Syn);               % original backbone configs

tic

D = NaN * ones(3,6);

switch length(Move),                            % how many sugars need moving

case 1,

  for f1 = 1:length(File),
    for i = 1:length(File(f1).NT),
      if File(f1).NT(i).Base == Motif.NT(Move(1)).Base,
        A = File(f1).NT(i);                     % base from 3D structure
        Sugar = A.Sugar(1:12,:) - ones(12,1)*A.Center;  % subtract center
        Sugar = Sugar * A.Rot;                  % rotate appropriately

        N1 = Motif.NT(Move(1));
        F.NT(Move(1)).Sugar = (ones(12,1)*N1.Center + Sugar*N1.Rot');
        F.NT(Move(1)).Syn   = File(f1).NT(i).Syn;

        D = NaN * ones(3,6);

        for m = Move,
          for b = 1:3,
            if b == m,                          % base the sugar is attached to
              D(m,b) = Inf;                     % any distance is OK
            else
              D(m,2*b-1) = min(min(zDistance(F.NT(m).Sugar,F.NT(b).Sugar)));
              D(m,2*b) = min(min(zDistance(F.NT(m).Sugar,F.NT(b).Fit)));
            end
          end
        end

        att = att + 1;

        us = unique(Sugar(1:12,:),'rows');           % some sugars use N1/N9

if max(max(Sugar)) > 10,
  [f1 i]
  Sugar
end

        Same = all(OrigConfig == cat(1,F.NT.Syn));

        if min(min(D)) > max(Mod(:,Same+5)) && length(us(:,1)) == 12,   
                                                % better min distance

figure(4)
clf
zDisplayNT(F,[2 3 1],VP);
text(0,2,0,num2str(min(min(D))))
axis off
view(2)
drawnow

          Motif = F;                      % store this one

          c = c + 1;
          Mot{c} = F;
          Mod(c,Same+5) = min(min(D));
          Mod(c,1) = f1;
          Mod(c,3) = i;
        end
      end
if min(min(D)) > BigDistCutoff && Same == 1,
  break
end

    end
if min(min(D)) > BigDistCutoff && Same == 1,
  break
end
  end

case 2,

  for f1 = 1:length(File),
    for i = 1:length(File(f1).NT),
      if File(f1).NT(i).Base == Motif.NT(Move(1)).Base,
        A = File(f1).NT(i);                     % base from 3D structure
        Sugar = A.Sugar(1:12,:) - ones(12,1)*A.Center;  % subtract center
        Sugar = Sugar * A.Rot;                  % rotate appropriately

        us = unique(Sugar,'rows');      % some sugars use N1/N9

        N1 = Motif.NT(Move(1));
        F.NT(Move(1)).Sugar = (ones(12,1)*N1.Center + Sugar*N1.Rot');
        F.NT(Move(1)).Syn   = File(f1).NT(i).Syn;

        DD = NaN * ones(1,3);
        m = Move(1);                            % backbone being moved now
          for b = 1:3,
            if b == m,                          % base the sugar is attached to
              DD(b)   = Inf;
            else
              DD(b)   = min(min(zDistance(F.NT(m).Sugar,F.NT(b).Fit)));
            end
          end

        if length(us(:,1)) == 12 && min(DD) > 2,  
                                                 % all atoms are distinct
                                                 % bb is far from other bases
        for f2 = 1:length(File),
          for j = 1:length(File(f2).NT),
            if File(f2).NT(j).Base == Motif.NT(Move(2)).Base,

              A = File(f2).NT(j);                     % base from 3D structure
              Sugar = A.Sugar(1:12,:) - ones(12,1)*A.Center;  % subtract center
              Sugar = Sugar * A.Rot;                  % rotate appropriately

              N2 = Motif.NT(Move(2));
              F.NT(Move(2)).Sugar = (ones(12,1)*N2.Center + Sugar*N2.Rot');
              F.NT(Move(2)).Syn   = File(f2).NT(j).Syn;

if mod(c,10000) == 0,
  toc
  [c f1 f2 length(File)]
  tic
end

              D = NaN * ones(3,6);
              for m = Move,
                for b = 1:3,
                  if b == m,                   % base the sugar is attached to
                    D(m,2*b-1) = Inf;          % any distance is OK
                    D(m,2*b)   = Inf;
                  else
                    D(m,2*b-1) = min(min(zDistance(F.NT(m).Sugar,F.NT(b).Sugar))); 
                    D(m,2*b)   = min(min(zDistance(F.NT(m).Sugar,F.NT(b).Fit)));
                  end
                end
              end

              att = att + 1;

              if mod(att,100000) == 0,
                fprintf('Made %d attempts so far, up to %7.4f\n', att, max(Mod(:,6)));
                drawnow
              end

              us = unique(Sugar(1:12,:),'rows');     % some sugars use N1/N9

              Same = all(OrigConfig == cat(1,F.NT.Syn));

              if min(min(D)) > max(Mod(:,Same+5)) && length(us(:,1)) == 12,   
                                                    % better min distance
figure(4)
clf
zDisplayNT(F,[2 3 1],VP);
text(0,2,0,num2str(min(min(D))))
axis off
view(2)
drawnow

                Motif = F;                      % store this one

                c = c + 1;

                Mot{c} = F;
                Mod(c,Same+5) = min(min(D));
                Mod(c,1) = f1;
                Mod(c,2) = f2;
                Mod(c,3) = i;
                Mod(c,4) = j;

              end
            end
if (min(min(D)) > BigDistCutoff && Same == 1) || att >= MaxTry,
  break
end
          end
if (min(min(D)) > BigDistCutoff && Same == 1) || att >= MaxTry,
  break
end
        end
        end                                % if sugar has 12 atoms
      end
if (min(min(D)) > BigDistCutoff && Same == 1) || att >= MaxTry,
  break
end
    end
if (min(min(D)) > BigDistCutoff && Same == 1) || att >= MaxTry,
  break
end
  end

case 3,
  fprintf('Okay, this motif needs 3 sugars moved!\n');

  [S1B3 S2B3 S3B1 S1S3 S2S3 S2B1 S2S1]

  [y,i] = sort([S1B3 S2B3 S3B1 S1S3 S2S3 S2B1 S2S1]);

  if S1S3 < MinDistCutoff,
    Motif = zBackboneSubstitution(Motif,File,AltSugar,[1 3]);
    Move = setdiff(Move,[1 3]);
  end

  if S2S3 < MinDistCutoff,
    Motif = zBackboneSubstitution(Motif,File,AltSugar,[2 3]);
    Move = setdiff(Move,[2 3]);
  end

  if S2S1 < MinDistCutoff,
    Motif = zBackboneSubstitution(Motif,File,AltSugar,[2 1]);
    Move = setdiff(Move,[2 1]);
  end

  if S1B3 < MinDistCutoff,
    Motif = zBackboneSubstitution(Motif,File,AltSugar,[1]);
    Move = setdiff(Move,[1]);
  end

  if S2B3 < MinDistCutoff,
    Motif = zBackboneSubstitution(Motif,File,AltSugar,[2]);
    Move = setdiff(Move,[2]);
  end

  if S2B1 < MinDistCutoff,
    Motif = zBackboneSubstitution(Motif,File,AltSugar,[2]);
    Move = setdiff(Move,[2]);
  end

  if S3B1 < MinDistCutoff,
    Motif = zBackboneSubstitution(Motif,File,AltSugar,[3]);
    Move = setdiff(Move,[3]);
  end

  Move = [];

end

% ----------------------------------------------------- Display improvements

if Verbose > 2,

D = cat(1,Mod(:,1));

[y,k] = sort(-D);

for m = 1:length(k),
        figure(2)
        clf
        c = k(m);

        f1 =      Mod(c,1);
        f2 =      Mod(c,2);
        i  =      Mod(c,3);
        j  =      Mod(c,4);

        F = Motif;

        A = File(f1).NT(i);                     % base from 3D structure
        Sugar = A.Sugar - ones(13,1)*A.Center;  % subtract center
        Sugar = Sugar * A.Rot;                  % rotate appropriately

        N1 = Motif.NT(Move(1));
        F.NT(Move(1)).Sugar = (ones(13,1)*N1.Center + Sugar*N1.Rot');

        A = File(f2).NT(j);                     % base from 3D structure
        Sugar = A.Sugar - ones(13,1)*A.Center;  % subtract center
        Sugar = Sugar * A.Rot;                  % rotate appropriately

        N2 = Motif.NT(Move(2));
        F.NT(Move(2)).Sugar = (ones(13,1)*N2.Center + Sugar*N2.Rot');

        FN = File(f1).Filename;
        Ba = [File(f1).NT(i).Base File(f1).NT(i).Number];

        zDisplayNT(F,[2 3 1],VP);
        T = [Motif.Filename ' Model ' num2str(c) ' with sugar substituted from ' FN ' ' Ba];
        T = strrep(T,'_','-');
        title(T);
        zLinkFigures([1 2]);
        view(2)
        zWritePDB(F,[Motif.Filename '_Model_' num2str(c) '.pdb']);

end
end

else


end

% ------------------------------------------ Decide which motif to use

if length(Move) > 0,
  [z,d] = max(Mod(:,5));                   % best with different config
  [y,c] = max(Mod(:,6));                   % best with same config

[z d y c]

  if y >= 3 || y >= z - 0.5,
    Motif = Mot{c};
    m = y;
  else
    Motif = Mot{d};
    m = z;
    fprintf('*********************** Used a backbone conformation change\n');
  end

  fprintf('Improved minimum distance from %7.4f to %7.4f in %d attempts\n', InitDist, m, att);
end
