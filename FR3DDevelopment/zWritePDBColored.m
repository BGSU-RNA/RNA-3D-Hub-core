% zWritePDBColored writes the nucleotide information in File to a PDB file that Swiss PDB can read.  Nucleotides are colored according to Method and the file is given the name Filename, if specified.  Method 2 uses the vector S, coloring from red for the minimum value of S to blue for the maximum value of S.
% Examples:
% zWritePDBColored('354D')
% zWritePDBColored(File,1,'354D_colored.pdb');
% zWritePDBColored(File,2,[File.Filename '_colored_by_S.pdb'],S)

function [void] = zWritePDBColored(File,Method,Filename,S,CC)

if strcmp(class(File),'char'),
  Filename = File;
  File = zAddNTData(Filename);
end

if nargin < 2,
  Method  = 1;
end

N = length(File.NT);                           % number of nucleotides

colormap('default');
map = colormap;

switch Method
case 1,                                        % color by index
  FN = [File.Filename '_colored_by_index.pdb'];
  C = map(7+ceil(48*(1:N)/N),:);
case 2,                                        % color by user-supplied vector
  FN = [File.Filename '_colored.pdb'];
  mins = min(S);
  maxs = max(S);
  C = map(7+ceil(48*(S-mins)/(maxs-mins)),:);
end

if nargin == 5,
  C = CC;
end


if nargin < 3,
  Filename = FN;
end

% ---------------------------------------------- Write the PDB file

fid = fopen(Filename,'w');       % open for writing

a = 1;                                         % atom number

for n = 1:length(File.NT),
  a = zWriteNucleotidePDB(fid,File.NT(n),a);
end

% ---------------------------------------------- Write color information

fprintf(fid,'SPDBVi 1 1 1 0 1 0 1 1 0 1 0 1 1 0 0\n');

L = N;

  while L > 0,
    if L >= 20,
      fprintf(fid,'SPDBVf 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19\n');
      L = L - 20;
    else
      fprintf(fid,'SPDBVf ');
      for i = 1:L,
        fprintf(fid,'19 ');
      end
      fprintf(fid,'\n');
      L = 0;
    end
  end    

  % Write colors for each element

L = N;                                      % number of nucleotides
  a = 1;

  while L > 0,
    if L >= 4,
      fprintf(fid, 'SPDBVc ');
      fprintf(fid, '%4.2f ', C(a,:));
      fprintf(fid, '%4.2f ', C(a+1,:));
      fprintf(fid, '%4.2f ', C(a+2,:));
      fprintf(fid, '%4.2f ', C(a+3,:));
      fprintf(fid,'\n');
      L = L - 4;
      a = a + 4;
    else
      fprintf(fid, 'SPDBVc ');
      for i = 1:L,
        fprintf(fid,'%4.2f ', C(a,:));
        L = L - 1;
        a = a + 1;
      end
    end
  end

fprintf(fid,'\n');

fclose(fid)
