% zSymmetryUnits reads .pdb1, pdb2, etc. files to find additional interactions

% First issue:  Identify cases in which the .pdb file contains multiple
% biological units, so that the pdb1, pdb2, files contain strictly less
% data.  Only apply this procedure in cases in which it adds something.
% But how to tell?  Maybe John Westbrook can help.

FN = '2Y9H.pdb';        % ??? different files seem to have different structures
                        % asymmetric unit contains 8 biological units
N = 8;


FN = '1HR2.pdb';        % asymmetric unit contains two biological units,
                        % so just use it, don't use the others
N = 2;

clear File

for v = 1:N,
  Filename = [FN num2str(v)];
  F = zAddNTData(Filename);

  if v == 1,
    File.NT = F.NT;
    File.Edge = F.Edge;
    File.Filename = F.Filename;
  else
    File.NT = [File.NT F.NT];
    [s,t] = size(File.Edge);
    [ss,tt] = size(F.Edge);

%    File.Edge = [File.Edge zeros(s,tt); zeros(ss,t) F.Edge];

  end

end

FF.Filename = [File.Filename '_symmetry'];
FF.NT       = File.NT;
FF.NumNT    = length(File.NT);

zSaveNTData(FF);

FF = zGetNTData(FF.Filename,1,1);
FF = zBackboneContinuity(FF);

figure(1)
clf
VP.Sugar = 1;
VP.LabelBases = 0;
zDisplayNT(FF,'all',VP);

figure(2)
clf
zCircularDiagram(FF);


break

  a = ['http://www.rcsb.org/pdb/files/' Filename];

    c = urlread(a);
    L = length(c);
    T = reshape(c',81,[])';
    T = T(:,1:80);

    fid = fopen(['PDBFiles' filesep Filename],'w');
    for i = 1:length(T(:,1)),
      fprintf(fid,'%s\n',T(i,:));
    end
    fclose(fid);
