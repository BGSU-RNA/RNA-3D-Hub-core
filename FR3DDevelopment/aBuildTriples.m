% Creates models for all possible triples and saves them as pdb and png
% files. Also writes out a csv table with a summary.
function [] = aBuildTriples(Location)

destination = [Location filesep 'Models'];
if ~exist(destination,'dir')
    mkdir(destination);
end
tables = [Location filesep 'Tables'];
if ~exist(tables,'dir')
    mkdir(tables);
end
fid = fopen([tables filesep 'Triple_models.csv'],'w');

Letter = 'ACGU';
% cWW tWW cHW tHW cSW tSW cWH tWH cHH tHH cSH tSH - correspond to the rows
first = [1 2 -3 -4 -5 -6 3 4 7 8 -9 -10];
% cHW tHW cHH tHH cHS tHS cSW tSW cSH tSH cSS tSS - correspond to the columns
second = [-3 -4 -7 8 9 10 -5 -6 -9 -10 11 12];
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
                            [F,message] = aMakeTriple(BP1,BP2,a,b,c);
                            if ~isempty(F)
                                created = created + 1;
                                clf ('reset');
                                VP.Sugar = 1;
                                VP.AtOrigin = 1;
                                zDisplayNT(F,[2 1 3],VP)
                                Title = 'Model ';
                                Title = [Title F.NT(1).Base F.NT(1).Number '-'];
                                Title = [Title F.NT(2).Base F.NT(2).Number '-'];
                                Title = [Title F.NT(3).Base F.NT(3).Number ' '];
                                Title = [Title strrep(BP1,' ','') '-' strrep(BP2,' ','')];
                                title(Title);
                                axis off
                                view(2)
                                Filename = [destination filesep name '_' Letter(a) Letter(b) Letter(c) '.png'];
                                saveas(gcf,Filename,'png');
                                Filename = [destination filesep name '_' Letter(a) Letter(b) Letter(c) '.pdb'];
                                zWritePDB(F,Filename,F.NT(2).Rot,F.NT(2).Fit(1,:));
                                fprintf(fid,'%s_%s,%s%s%s,Yes\n',BP1,BP2,Letter(a),Letter(b),Letter(c));
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
end


function [File,message] = aMakeTriple(Pair1,Pair2,Base1,Base2,Base3)

    clear File
    message = '';
    Letters = 'ACGU';
    % retrieve examplar basepairs
    [N1,N2] = zGetExemplar(Pair1,Base1,Base2);
    [M2,M3] = zGetExemplar(Pair2,Base2,Base3);
    
    if ~isempty(N1.Code) && ~isempty(N2.Code) && ~isempty(M2.Code) && ~isempty(M3.Code),

        if 0 > 1,
          figure(2)
          clf
          VP.Sugar = 1;
          VP.AtOrigin = 1;
          zPlotOneNT(M2,VP)
          hold on
          zPlotOneNT(M3,VP)
          view(2)
          figure(1)
        end

        Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
        Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

        % fprintf('%s%s%s %5s %5s\n',N1.Base, N2.Base, M3.Base, zEdgeText(

        % rotate the second pair so that the middle bases (M2, N2) superimpose

        [s,t] = size(M2.Fit);

        M2.Fit = (M2.Fit - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;

        [s,t] = size(M3.Fit);

        M3.Fit = (M3.Fit - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;

        [s,t] = size(M2.Sugar);
        M2.Sugar = (M2.Sugar - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;
        M3.Sugar = (M3.Sugar - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;

        % M3.Rot = M3.Rot * M2.Rot' * N2.Rot;
        M3.Rot = N2.Rot * M2.Rot' * M3.Rot;

        % add the three nucleotides to a File data structure
        File.NT(1) = N1;
        File.NT(1).Chain = 'A';
        File.NT(2) = N2;
        File.NT(1).Chain = 'B';
        File.NT(3) = M3;
        File.NT(1).Chain = 'C';        
        File.NumNT = 3;

        for i = 1:3,
            File.NT(i).Center = mean(File.NT(i).Fit(1:Lim(1,File.NT(i).Code),:));
        end
        
    else
        File = [];
        if isempty(N1.Code) || isempty(N2.Code)
            message = sprintf('%s %s%s does not exist',Pair1,Letters(Base1),Letters(Base2));
        end        
        if isempty(M2.Code) || isempty(M3.Code)
            if ~isequal(message,'')
                message = [message ' and '];
            end
            message = [message sprintf('%s %s%s does not exist',Pair2,Letters(Base2),Letters(Base3))];
        end
    end

end
    

% function [File,message] = zMakeTriple(Pair1,Pair2,Base1,Base2,Base3)
% 
% clear File
% message = '';
% 
% % retrieve examplar basepairs
% 
% [N1,N2] = zGetExemplar(Pair1,Base1,Base2);
% [M2,M3] = zGetExemplar(Pair2,Base2,Base3);
% 
% if ~isempty(N1.Code) && ~isempty(N2.Code) && ~isempty(M2.Code) && ~isempty(M3.Code),
% 
%     if 0 > 1,
%         figure(2)
%         clf
%         VP.Sugar = 1;
%         VP.AtOrigin = 1;
%         zPlotOneNT(M2,VP)
%         hold on
%         zPlotOneNT(M3,VP)
%         view(2)
%         figure(1)
%     end
% 
%     Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
%     Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen
% 
%     % fprintf('%s%s%s %5s %5s\n',N1.Base, N2.Base, M3.Base, zEdgeText(
% 
%     % rotate the second pair so that the middle bases (M2, N2) superimpose
% 
%     [s,t] = size(M2.Fit);
% 
%     M2.Fit = (M2.Fit - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;
% 
%     [s,t] = size(M3.Fit);
% 
%     M3.Fit = (M3.Fit - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;
% 
%     [s,t] = size(M2.Sugar);
%     M2.Sugar = (M2.Sugar - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;
%     M3.Sugar = (M3.Sugar - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;
% 
%     % M3.Rot = M3.Rot * M2.Rot' * N2.Rot;
%     M3.Rot = N2.Rot * M2.Rot' * M3.Rot;
% 
%     % add the three nucleotides to a File data structure
% 
%     File.NT(1) = N1;
%     File.NT(1).Chain = 'A';
%     File.NT(2) = N2;
%     File.NT(2).Chain = 'B';
%     %M3.Center =
%     File.NT(3) = M3;
%     File.NT(3).Chain = 'C';
%     %File.NT(4) = M2;
%     File.NumNT = 3;
% 
%     for i = 1:3,
%         File.NT(i).Center = mean(File.NT(i).Fit(1:Lim(1,File.NT(i).Code),:));
%     end
% 
% else
%     File = [];
% end
% 
% end