% xAnnotate(File) annotates a 3D structure with secondary structure elements

File = zAddNTData('2AVY');

File = xAnnotateWithKnownMotifs(File,1,0,{'Pseudoknot','Helix-----','Hairpin-----','Internal_loop','Bulge_loop','Junction_loop'});

clear T

zInteractionsExcel;

C.Name = '5 prime end';
C.Index = 1;
C.Number = 'x';

for i = 1:length(File.NT),
  if length(File.Nucl(i).Motif) > 0,
    C = File.Nucl(i).Motif(1);
  elseif isempty(strfind(C.Name,'Helix')) && isempty(strfind(C.Name,'Pseudo')),
    File.Nucl(i).Motif(1) = C;
  end
end  

C.Name = '3 prime end';
C.Index = 1;
C.Number = 'x';

for i = length(File.NT):-1:1,
  if length(File.Nucl(i).Motif) > 0,
    C = File.Nucl(i).Motif(1);
  elseif isempty(strfind(C.Name,'Helix')) && isempty(strfind(C.Name,'Pseudo')),
    File.Nucl(i).Motif(1) = C;
  end
end  

for i = 1:length(File.NT),
  for j = 1:length(File.Nucl(i).Motif),
    T{i+1,14+j} = strrep(File.Nucl(i).Motif(j).Name,'-----','');
  end
end


T(:,[1 2 15:end])

%T(:,[1 2 15])


