
function [File1,i1,File2,i2] = zReadJesseAlignment(Molecule)

switch Molecule

case '16S',  % ----------------------------------------------- 16S alignment

  [n,t,r] = xlsread('Stombaugh_et_al_Sup_Mat_S9.xls','16S_3D_Struct_Aln');

  cols = [5 7 18 20];
  startrow = 3;
  Chain1 = 'A';
  Chain2 = 'A';

  File1 = zAddNTData('2avy');
  File2 = zAddNTData('1j5e');

case '5S',  % ----------------------------------------------- 5S alignment

  [n,t,r] = xlsread('Stombaugh_et_al_Sup_Mat_S9.xls','5S_3D_Struct_Aln');

  cols = [5 7 18 20];
  startrow = 3;
  Chain1 = 'A';
  Chain2 = 'B';

  File1 = zAddNTData('2aw4');
  File2 = zAddNTData('2j01');

case '23S',  % ----------------------------------------------- 23S alignment

  [n,t,r] = xlsread('Stombaugh_et_al_Sup_Mat_S9.xls','23S_3D_Struct_Aln');

  cols = [5 7 18 20];
  startrow = 3;
  Chain1 = 'B';
  Chain2 = 'A';

  File1 = zAddNTData('2aw4');
  File2 = zAddNTData('2j01');

end


[i1,i2] = LookUpIndices(File1,File2,Chain1,Chain2,r,cols,startrow);

if 0 > 1,
for m = 1:20,
  fprintf('Jesse aligns %1s%4s with %1s%4s\n', File1.NT(i1(m)).Base, File1.NT(i1(m)).Number, File2.NT(i2(m)).Base, File2.NT(i2(m)).Number);
end
end

% ---------------------------------------------- LookUpIndices routine

function [i1,i2] = LookUpIndices(File1,File2,Chain1,Chain2,r,cols,startrow);

indices = [];

for a = startrow:length(r(:,1)),

if 0 > 1,
fprintf('Row %d\n', a);
r{a,cols(1)}
r{a,cols(2)}
r{a,cols(3)}
r{a,cols(4)}
end

  r1 = MakeString(r{a,cols(1)});
  r2 = MakeString(r{a,cols(2)});
  r3 = MakeString(r{a,cols(3)});
  r4 = MakeString(r{a,cols(4)});

  if ~isempty(r1) && ~isempty(r3),            % row a has an alignment btw 1,3
    i = zIndexLookup(File1,r1,Chain1);
    k = zIndexLookup(File2,r3,Chain2);
    indices = [indices; [i k]];
  end

  if ~isempty(r2) && ~isempty(r4),            % row a has an alignment btw 2,4
    j = zIndexLookup(File1,r2,Chain1);
    m = zIndexLookup(File2,r4,Chain2);
    indices = [indices; [j m]];               % append these matches
  end
end

indices = unique(indices,'rows');             % remove redundant aligns

i1 = indices(:,1);
i2 = indices(:,2);

% ---------------------------------------------- Get a string routine

function [s] = MakeString(r)

if isnan(r),
  s = '';
elseif strcmp(class(r),'double'),
  s = num2str(r);
elseif strcmp(class(r),'char'),
  s = r;
end
