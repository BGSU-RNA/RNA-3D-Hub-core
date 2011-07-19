% zConvertSponerPDB makes small changes to the PDB files written by Jiri Sponer

% nam = 'average_var_1_prot_45-45.5ns'; zConvertSponerPDB(nam); clf; zDisplayNT([nam '_fixed']); zShowInteractionTable([nam '_fixed'],'all');

function [void] = zConvertSponerPDB(filename)

if isempty(strfind(lower(filename),'.pdb'))
  file     = filename;
  filename = [filename '.pdb'];
else
  file = strrep(filename,'.pdb','');
  file = strrep(file,'.PDB','');
end

inid = fopen(['PDBFiles' filesep filename],'r');
outid = fopen(['PDBFiles' filesep file '_fixed.pdb'],'w');

      L = 1;
      while L > -1
        L = fgets(inid);
        if L > -1
          L = strrep(L,'RA ','  A');
          L = strrep(L,'RC ','  C');
          L = strrep(L,'RG ','  G');
          L = strrep(L,'RU ','  U');
          L = strrep(L,'RA3','  A');
          L = strrep(L,'RC3','  C');
          L = strrep(L,'RG3','  G');
          L = strrep(L,'RU3','  U');
          L = strrep(L,'RA5','  A');
          L = strrep(L,'RC5','  C');
          L = strrep(L,'RG5','  G');
          L = strrep(L,'RU5','  U');
          L = strrep(L,'RAP','  A');
          L = strrep(L,'RCP','  C');
          L = strrep(L,'RGP','  G');
          L = strrep(L,'RUP','  U');

          fprintf(outid,'%s',L);
        end
      end

fclose(inid);
fclose(outid);
