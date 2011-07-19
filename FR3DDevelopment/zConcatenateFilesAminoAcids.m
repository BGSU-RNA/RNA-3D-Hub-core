
function [F] = zConcatenateFilesAminoAcids(F1,F2)

% if File is a text string (filename), load the file

if strcmp(class(F1),'char'),
  Filename = F1;
  F1 = zGetAAData(Filename,0);
end

if strcmp(class(F2),'char'),
  Filename = F2;
  F2 = zGetAAData(Filename,0);
end

F.Filename = [F1.Filename '_' F2.Filename '_AA'];

U1 = upper(unique(cat(2,F1.AA.Chain)));
U2 = upper(unique(cat(2,F2.AA.Chain)));
U  = U1;                                  % accumulate used chains

Alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';      % unused chain identifiers

for i = 1:length(U1),
  j = find(U1(i) == Alph);
  if ~isempty(j),
    Alph = Alph([1:(j(1)-1) (j(1)+1):end]);
  end
end

for u = length(U2):-1:1,
  if any(U2(u) == U),
    c = Alph(1);                          % new chain letter to use
    Alph = Alph(2:end);                   % remove this letter
    U    = [U c];                         % add this to used list
    i = find(cat(1,F2.AA.Chain) == U2(u));
    for ii = 1:length(i),
      F2.AA(i(ii)).Chain = c;
    end
  else
    j = find(U2(u) == Alph);
    if ~isempty(j),
      Alph = Alph([1:(j(1)-1) (j(1)+1):end]);
    end
    U    = [U U2(u)];                     % add this to used list
  end
end

F.AA = [F1.AA F2.AA];
F.NumAA = length(F.AA);

%F.Info.Resolution = max(F1.Info.Resolution,F2.Info.Resolution);
%F.Info.Descriptor = [F1.Info.Descriptor '_' F2.Info.Descriptor];
%F.Info.ExpTechnique = [F1.Info.ExpTechnique '_' F2.Info.ExpTechnique];
%F.Info.ReleaseDate = [F1.Info.ReleaseDate '_' F2.Info.ReleaseDate];
%F.Info.Author = [F1.Info.Author '_' F2.Info.Author];
%F.Info.Keywords = [F1.Info.Keywords '_' F2.Info.Keywords];
%F.Info.Source = [F1.Info.Source '_' F2.Info.Source];

%zSaveNTData(F)

%F = zAddNTData(F.Filename,0,[],1);


