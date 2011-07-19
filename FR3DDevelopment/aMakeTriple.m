
% ======================================================================

function [File,message,AltSugar] = aMakeTriple(Pair1,Pair2,Base1,Base2,Base3)

    load PairExemplars_for_triple_models

    clear File
    message = '';
    Letters = 'ACGU';
    % retrieve examplar basepairs
    [N1,N2] = zGetExemplar(Pair1,Base1,Base2,Exemplar);
    [M2,M3] = zGetExemplar(Pair2,Base2,Base3,Exemplar);
    
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

% No need to create nucleotide M2 since N2 is used instead

%        M2.Fit = (M2.Fit - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;

        [s,t] = size(M3.Fit);

        M3.Fit = (M3.Fit - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;

        [s,t] = size(M2.Sugar);
        AltSugar = (M2.Sugar - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;
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
        AltSugar = [];
    end

end
