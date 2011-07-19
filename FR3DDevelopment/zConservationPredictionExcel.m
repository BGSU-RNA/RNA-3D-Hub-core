
clear Te
c = 1;

for r = 1:length(Data(:,1)),
  f = Data(r,1);
  
  if f > 0,

    i = Data(r,2);
    E = fix(File(f).Edge(i,:));
    B = fix(File(f).BasePhosphate(i,:));
    B(i) = 0;                              % ignore self interactions

    Te{c,1} = c;
    Te{c,2} = File(f).Filename;
    Te{c,3} = [File(f).NT(i).Base File(f).NT(i).Number];

    if ~isempty(File(f).Nucl(i).Motif),
      Te{c,4} = strrep(File(f).Nucl(i).Motif(1).Name,'-','');
      Te{c,4} = strrep(Te{c,4},'_',' ');
    else
      Te{c,4} = '';
    end

    j = find(abs(E) == 1);    
    if ~isempty(j),
      newt = [File(f).NT(j).Base File(f).NT(j).Number];
    else
      newt = '';
    end
    Te{c,5} = newt;

    j = find((abs(E) > 1) .* (abs(E) < 13));    
    newt = '';
    if ~isempty(j),
      for jj = 1:length(j),
        newt = [newt File(f).NT(j(jj)).Base File(f).NT(j(jj)).Number ' ' zEdgeText(E(j(jj))) ' '];
      end
    end
    Te{c,6} = newt;

    j = find((abs(B) > 0) .* (abs(B) < 30));    
    newt = '';
    if ~isempty(j),
      for jj = 1:length(j),
        newt = [newt File(f).NT(j(jj)).Base File(f).NT(j(jj)).Number zBasePhosphateText(B(j(jj))) ' '];
      end
    end
    Te{c,7} = newt;
    
    Te{c,8} = Data(r,7);
    Te{c,9} = 100*Data(r,18);
    Te{c,10} = Data(r,3);
    Te{c,11} = Data(r,19);
    c = c + 1;
  end
end

H = {'Index','File','Nucleotide','Secondary structure','cWW pair','non-cWW pairs','Base-phosphate','Near protein','Stacking partner conservation','Conservation','Predicted Conservation'};

Te = [H; Te];

xlswrite('Conservation_and_interaction_data_E_coli.xls',Te);
