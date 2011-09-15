
loopType = 'HL';                            % hairpins only
loopType = 'IL';                           % internal loops only

ModelFNs = dir(['Models' filesep]);
keep = zeros(length(ModelFNs));
for m = 1:length(ModelFNs),
   if (length(ModelFNs(m).name) > 2),
      if (strcmp(ModelFNs(m).name(1:2),loopType) && ~strcmp(ModelFNs(m).name(4),'M') && ~strcmp(ModelFNs(m).name(4),'.')), 
          keep(m) = 1;
          ModelFNs(m).name = strrep(ModelFNs(m).name,'.txt','');
      end
   end 
end
ModelFNs = ModelFNs(find(keep));
JAR3D_path;
for m = 1:length(ModelFNs),
  MN = ModelFNs(m).name;
  MFN = [MN '.txt'];
  SFN = ['RS_' MN(1:6) '.fasta'];
  Scores = JAR3DMatlab.MotifParseSingle(pwd,SFN,MFN);
  OFN = ['Models' filesep 'Emp. Distributions' filesep MN(1:6) '.txt'];
%  save(OFN,'Scores','-mat')
  fid = fopen(OFN,'w');
  for i = 1:length(Scores)
    fprintf(fid, '%f\n', Scores(i));
  end
  fclose(fid);
end
