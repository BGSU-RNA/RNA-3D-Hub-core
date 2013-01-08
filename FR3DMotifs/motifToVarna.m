function [Text] = motifToVarna(Search, Location, execute)

    javaCommand = 'java -Djava.awt.headless=true -cp ../VARNAv3-7.jar fr.orsay.lri.varna.applications.VARNAcmd';
    destination = [Location filesep '2ds'];
    if ~exist(destination,'dir'), mkdir(destination); end

    if nargin<3
        execute = 1; % run unix commands
    end

    Text = [];

    nts = get_sequence_consensus(Search);
    str(1:length(nts)-2) = '.';

    Text = [Text sprintf('-sequenceDBN "%s" ',nts)];
    Text = [Text sprintf('-structureDBN "(%s)" ',str)];
    Text = [Text '-baseNum "#ffffff" '];
    Text = [Text sprintf('-auxBPs "')];

    % (1,13):edge5=wc,edge3=wc,stericity=cis;

    if isfield(Search,'consensusEdge')
        Edge = Search.consensusEdge;
    else
        Edge = pConsensusInteractions(Search);
    end

    EdgeAbsFix = fix(abs(Edge));

    N = length(Edge);

    for i = 1:N
        for j = (i+1):N

            if EdgeAbsFix(i,j) <= 12 && EdgeAbsFix(i,j) > 0

                if mod(EdgeAbsFix(i,j),2) ~= 0
                    stericity = 'cis';
                else
                    stericity = 'trans';
                end

                int = zEdgeText(Edge(i,j));
                edge5 = lower(int(2));
                edge3 = lower(int(3));
                if edge5 == 'w'
                    edge5 = 'wc';
                end
                if edge3 == 'w'
                    edge3 = 'wc';
                end
                Text = [Text sprintf('(%i,%i):edge5=%s,edge3=%s,stericity=%s;',i,j,edge5,edge3,stericity)];

            end
        end
    end
    Text = [Text '"'];

    Text = sprintf('%s %s -o %s/%s.png',javaCommand,Text,destination,Search.FileName(1:end-4));

    if execute == 1
        unix(Text);
    else
        fprintf('%s\n',Text);
    end

end

function [consensus] = get_sequence_consensus(Search)

    [r,c] = size(Search.Candidates);
    ind = c-1;

    nts = zeros(r,ind);

    for i = 1:r

        if isfield(Search,'QIndex')
            targetFile = Search.File(Search.QIndex((Search.Candidates(i,end))));
        else
            targetFile = Search.File(Search.Candidates(i,end));
        end

        indices = Search.Candidates(i,1:end-1);
        nts(i,1:ind) = [targetFile.NT(indices).Base];
    end

    nts = char(nts);

    consensus = '';
    for i = 1:ind
        consensus(i) = analyze_sequence_column(nts(:,i));
    end

end

function [consensus] = analyze_sequence_column(nts)

    nts = reshape(nts,1,[]);

    if isempty(strfind(nts,'A'))
        A = 0;
    else
        A = 1;
    end
    if isempty(strfind(nts,'C'))
        C = 0;
    else
        C = 1;
    end
    if isempty(strfind(nts,'G'))
        G = 0;
    else
        G = 1;
    end
    if isempty(strfind(nts,'U'))
        U = 0;
    else
        U = 1;
    end

    if A == 0 && G == 0 && C == 1 && U == 1
        consensus = 'Y'; % pyrimidines
    elseif A == 1 && G == 1 && C == 0 && U == 0
        consensus = 'R'; % purines
    elseif A == 1 && G == 1 && C == 1 && U == 1
        consensus = 'N'; % any nucleotide
    elseif A == 1 && G == 0 && C == 0 && U == 0
        consensus = 'A';
    elseif A == 0 && G == 1 && C == 0 && U == 0
        consensus = 'G';
    elseif A == 0 && G == 0 && C == 1 && U == 0
        consensus = 'C';
    elseif A == 0 && G == 0 && C == 0 && U == 1
        consensus = 'U';
    elseif A == 1 && G == 0 && C == 0 && U == 1
        consensus = 'W'; % weak
    elseif A == 0 && G == 1 && C == 1 && U == 0
        consensus = 'S'; % strong
    elseif A == 0 && G == 1 && C == 0 && U == 1
        consensus = 'K'; % keto
    elseif A == 1 && G == 0 && C == 1 && U == 0
        consensus = 'M'; % amino
    elseif A == 1 && G == 1 && C == 0 && U == 1
        consensus = 'D'; % not C
    elseif A == 1 && G == 1 && C == 1 && U == 0
        consensus = 'V'; % not U
    elseif A == 1 && G == 0 && C == 1 && U == 1
        consensus = 'H'; % not G
    elseif A == 0 && G == 1 && C == 1 && U == 1
        consensus = 'B'; % not A
    else
        error('Unexpected case');
    end

end


%     if isfield(Search,'QIndex')
%         targetFile = Search.File(Search.QIndex((Search.Candidates(1,end))));
%     else
%         targetFile = Search.File(Search.Candidates(1,end));
%     end
%
%
%     indices = Search.Candidates(1,1:end-1);
%
%     nts = [targetFile.NT(indices).Base];
%     str(1:length(targetFile.NT(indices))-2) = '.';

    % Text = [Text sprintf('<param name="sequenceDBN%i"  value="%s" />\n',order,nts)];
    % Text = [Text sprintf('<param name="structureDBN%i"  value="(%s)" />\n',order,str)];
    % Text = [Text sprintf('<param name="title" value="LIBI%04d" />\n',order)];
    % Text = [Text sprintf('<param name="auxBPs%i" value="',order)];
    % Text = [Text '"/>\n'];
% <parameter name="applyBasesStyle5on"
% value="1,6,7,10"/>
% <parameter name="basesStyle3"
% value="fill=#FF0000,outline=#00FF00"/>
