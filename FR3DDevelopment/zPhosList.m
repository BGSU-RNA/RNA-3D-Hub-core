% zPhosList lists base-phosphate interaction information

% if Param(1) = 1, show near interactions.  0, don't
% if Param(2) = 1, show self interaction.  0, don't

function [void] = zPhosList(File,D,Param)

if nargin < 3,
  Param = [1 1];
end

D = sortrows(D,[1 2 3 5 7]);

L = {'A','C','G','U'};
pn  = {'O5*','O3*','O1 ','O2 '};              % names of phosphate oxygens

if Param(1) == 0,
  r = find(D(:,5) < 100);
  D = D(r,:);                            % remove near interactions
end

if Param(2) == 0,
  r = find(D(:,2) ~= D(:,3));
  D = D(r,:);
end

for g = 1:19,

  c = 1;
  r = find(D(:,23) == g);                  % rows for this interaction

  if length(r) > 0,

    DD = D(r,:);

    fprintf('Interaction list for %s %s interaction.  %d found.\n\n', L{DD(1,4)}, zBasePhosphateText(g),length(r));

    clear Tex
    clear Numb

    j = 1;
    while j <= length(r),
      Tex = sprintf('File %s  ', File(DD(j,1)).Filename);
      Tex = [Tex sprintf('Base %s%4d Phosphate %4d ', L{DD(j,4)}, DD(j,21), DD(j,22))];

      switch DD(j,4),
        case 1,                         % Base A
              hn  = {'H2','H8','1H6','2H6'}; % names of the base hydrogens
        case 2,                         % Base C
              hn  = {'H6','H5','1H4','2H4'}; % names of the base hydrogens
        case 3,                         % Base G
              hn  = {'H1','H8','1H2','2H2'}; % names of the base hydrogens
        case 4,                         % Base U
              hn  = {'H5','H3','H6'}; % names of the base hydrogens
      end
  
      Tex = [Tex sprintf('Donor %4s Oxygen %4s %5s ', hn{DD(j,6)}, pn{DD(j,7)}, zBasePhosphateText(DD(j,5)))];

      Numb(c,1) = DD(j,7);
      if DD(j,5) < 100,
        Numb(c,2) = 1;
      else
        Numb(c,2) = 0;
      end
      Numb(c,3) = 0;

      j = j + 1;

      done = 0;
      while done == 0 && j <= length(r),
        if all(DD(j-1,[1 2 3]) == DD(j,[1 2 3])),
          Tex = [Tex sprintf('| Donor %4s Oxygen %4s %5s ', hn{DD(j,6)}, pn{DD(j,7)}, zBasePhosphateText(DD(j,5)))];
          if DD(j,5) < 100, 
            Numb(c,2) = Numb(c,2) + 1;
            Numb(c,3) = DD(j,7);
          end
          j = j + 1;
        else
          done = 1;
        end
      end

      Text{c} = Tex;
      c = c + 1;
    end

    [y,i] = sortrows(Numb,[1 2 3]);
    Text = Text(i);
    for j = 1:length(Text),
      fprintf('%s\n', Text{j});
    end
    fprintf('\n')

  end
end

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
          %21  nucleotide number of base
          %22  nucleotide number of phosphate donor
          %23  (to be added below) classification of this interaction
