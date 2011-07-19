% xAnnotate annotates a 3D structure file with standard secondary structure elements.  The result is a new field File.Nucl whose first entry is the correct one

% File = zAddNTData('2avy');

File = xAnnotateWithKnownMotifs(File,1,0,{'Pseudoknot','Helix-----','Hairpin-----','Internal_loop','Bulge_loop','Junction_loop'});

clear T

for f = 1:length(File),

  %  zInteractionsExcel;

  C.Name = '5 prime end';
  C.Index = 1;
  C.Number = 'x';

  for i = 1:length(File(f).NT),
    if length(File(f).Nucl(i).Motif) > 0,
      C = File(f).Nucl(i).Motif(1);
    elseif isempty(strfind(C.Name,'Helix')) && isempty(strfind(C.Name,'Pseudo')),
      File(f).Nucl(i).Motif(1) = C;
    end
  end  

  C.Name = '3 prime end';
  C.Index = 1;
  C.Number = 'x';

  for i = length(File(f).NT):-1:1,
    if length(File(f).Nucl(i).Motif) > 0,
      C = File(f).Nucl(i).Motif(1);
    elseif isempty(strfind(C.Name,'Helix')) && isempty(strfind(C.Name,'Pseudo')),
      File(f).Nucl(i).Motif(1) = C;
    end
  end  

  for i = 1:length(File(f).NT),
    for j = 1:length(File(f).Nucl(i).Motif),
      T{i+1,14+j} = strrep(File(f).Nucl(i).Motif(j).Name,'-----','');
    end
  end

  %T(:,[1 2 15:end])

  %T(:,[1 2 15])
end

