% Creates models for all possible triples and saves them as pdb and png
% files. Also writes out a csv table with a summary.

%function [] = aBuildTriples(Location)

Location = 'Triples';                         % use current directory

if ~exist('File'),
  File = zAddNTData('Nonredundant_2A_2011-06-25_list');
end

destination = [Location filesep 'Models'];
if ~exist(destination,'dir')
    mkdir(destination);
end
tables = [Location filesep 'Tables'];
if ~exist(tables,'dir')
    mkdir(tables);
end

Letter = 'ACGU';
% cWW tWW cHW tHW cSW tSW cWH tWH cHH tHH cSH tSH - correspond to the rows
first = [1 2 -3 -4 -5 -6 3 4 7 8 -9 -10];
% cHW tHW cHH tHH cHS tHS cSW tSW cSH tSH cSS tSS - correspond to the columns
second = [-9 -10 11 12 -3 -4 -7 8 9 10 -5 -6];

switch processor,
case 1,
  first = [4 -6 7 8 -9 -10];
  diary TripleBackboneProc1b.txt
  fid = fopen([tables filesep 'Triple_models_proc1c.csv'],'w');
case 2,
  first = [-4 -5];
  diary TripleBackboneProc2b.txt
  fid = fopen([tables filesep 'Triple_models_proc2c.csv'],'w');
case 3,
  first = [1 2 -3 3];
  diary TripleBackboneProc3b.txt
  fid = fopen([tables filesep 'Triple_models_proc3c.csv'],'w');
case 4,
  first = [1 2 -3 -4 -5 -6 3 4 7 8 -9 -10];
  second = [-9 9];
  diary TripleBackboneProc4b.txt
  fid = fopen([tables filesep 'Triple_models_proc4b.csv'],'w');
otherwise,
  % cWW tWW cHW tHW cSW tSW cWH tWH cHH tHH cSH tSH - correspond to the rows
  first = [1 2 -3 -4 -5 -6 3 4 7 8 -9 -10];
  % cHW tHW cHH tHH cHS tHS cSW tSW cSH tSH cSS tSS - correspond to the columns
  second = [-9 -10 11 12 -3 -4 -7 8 9 10 -5 -6];
  fid = fopen([tables filesep 'Triple_models_proc.csv'],'w');
end

% M = zeros(12,12);

impossible = 0;
created    = 0;

for i = first
    BP1 = zEdgeText(i);
    BP1 = BP1(1:end-1);
    for j = second
        if ~(find(first==i)>6 && find(second==j)<7)
            BP2 = zEdgeText(j);
            BP2 = BP2(1:end-1);
            if ~isempty(intersect(i,[1 2 -3 -4 -5 -6])) && j == 12 % NB! cWW,tWW,cHW,tHW,cSW,tSW +_tsS
                BP2 = 'tsS';
            end                        
            name = ['Triple_' BP1 '_' BP2];
            if BP1(3) ~= BP2(2),
                fprintf('Making %s - %s triples\n', BP1, BP2);
                for a = 1:4,
                    for b = 1:4,
                        for c = 1:4,
                            %                             [F] = zMakeTriple(BP1,BP2,a,b,c);
                            [F,message,AltSugar] = aMakeTriple(BP1,BP2,a,b,c);

                            if ~isempty(F)

F.NT(1).Sugar = F.NT(1).Sugar(1:12,:);
F.NT(2).Sugar = F.NT(2).Sugar(1:12,:);
F.NT(3).Sugar = F.NT(3).Sugar(1:12,:);

                                Title = 'Model ';
                                Title = [Title F.NT(1).Base F.NT(1).Number '-'];
                                Title = [Title F.NT(2).Base F.NT(2).Number '-'];
                                Title = [Title F.NT(3).Base F.NT(3).Number ' '];
                                Title = [Title strrep(BP1,' ','') '-' strrep(BP2,' ','')];

                                figure(1)
                                clf ('reset');
                                VP.Sugar = 1;
                                VP.AtOrigin = 1;
                                VP.LabelBases = 0;
                                VP.LineThickness = 4; 
                                zDisplayNT(F,[2 1 3],VP);
                                title(Title);
                                axis off
                                view(2)
                                Filename = [destination filesep name '_' Letter(a) Letter(b) Letter(c) '_original.png'];
                                saveas(gcf,Filename,'png');


F.Filename = Title;
F = zBackboneSubstitution(F,File,AltSugar);



                                created = created + 1;
                                figure(1)
                                clf ('reset');
                                VP.Sugar = 1;
                                VP.AtOrigin = 1;
                                VP.LabelBases = 0;
                                zDisplayNT(F,[2 1 3],VP);
                                title(Title);
                                axis off
                                view(2)
                                Filename = [destination filesep name '_' Letter(a) Letter(b) Letter(c) '.png'];
                                saveas(gcf,Filename,'png');
                                Filename = [destination filesep name '_' Letter(a) Letter(b) Letter(c) '.pdb'];
                                zWritePDB(F,Filename,F.NT(2).Rot,F.NT(2).Fit(1,:));
                                fprintf(fid,'%s_%s,%s%s%s,Yes\n',BP1,BP2,Letter(a),Letter(b),Letter(c));

%disp('Pausing in aBuildTriplesNoClash');
%pause



                            else
                                impossible = impossible + 1;
                                fprintf(fid,'%s_%s,%s%s%s,%s\n',BP1,BP2,Letter(a),Letter(b),Letter(c),message);
                                disp(message);
                            end
                        end
                    end
                end
            end
        end
    end
end

fprintf('Created %i models, %i impossible\n',created, impossible);
close(gcf);

diary off
