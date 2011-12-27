% pMakeModelsFromLibrary is a script that:
% 1. Reads the .mat files in MotifLibrary
% 2. Finds coplanar and "loose" coplanar pairs (for the consensus basepairs)
% 3. Calls pMakeMotifModelFromSSF to make a model
% 4. Writes a .fasta file for each motif
% 5. Writes a .txt JAR3D file for each motif
% 6. Writes files of motif model names and fasta file names
% 7. Writes a file of signatures and reversed signatures

disp('Make sure the Matlab current folder has a MotifLibrary in it');

writeLD = 1;

group = cell(1);
sequence = cell(1);

if ~exist('loopType'),
    disp('Please specify a loop type, for example loopType = ''IL'';')
    break
end

JAR3D_path;

Param = [1 2 0 1 100 1];                   % See below

% Parameters stored in Param:
% Param(1) verbose
% Param(2) method to use for basepair isostericity
% Param(3) recognize extensible helices and model them as such
% Param(4) adjust substitution probabilities for long-range interactions
% Param(5) how far to look ahead for local basepair interactions
% Param(6) use near interactions

Prior = [.5 .5 .5 .5 0];              % Prior distribution for insertion bases

tableIndex = 0;                       % For use recording diagnostic information

% ----------------------------------------- Set variables

Verbose = Param(1);

% ----------------------------------------- Find .mat files for motifs

if ~(exist('MotifLibrary') == 7),        % if directory doesn't yet exist
    Filenames = dir;
    mkdir('MotifLibrary');
    for f = 1:length(Filenames),
        if ~isempty(strfind(Filenames(f).name,'.mat')),
            newname = strrep(Filenames(f).name,'Group',loopType);
            movefile(Filenames(f).name, ['MotifLibrary' filesep newname]);
        end
    end
end

% ----------------------------------------- Read file names from the library

Filenames = dir(['MotifLibrary' filesep]);

keep = [];                               % of all models, which to keep

for m = 1:length(Filenames),
    if (length(Filenames(m).name) > 2),
        if (Filenames(m).name(1:2) == loopType),
            keep(m) = 1;
            Filenames(m).name = strrep(Filenames(m).name,'.mat','');
        end
    end
end

Filenames = Filenames(find(keep));

% -----------------------------------------

CL = zClassLimits;                              % read ClassLimits matrix

if exist('PairExemplars.mat','file') > 0,
    load('PairExemplars','Exemplar');
else
    Exemplar = [];
end

% ----------------------------------------- Load each search, make a model

for m = 1:length(Filenames),
    MN = Filenames(m).name;
    FN = ['MotifLibrary' filesep MN '.mat'];
    load(FN,'Search','-mat')                             % Load search data

    %   if length(MN) < 9 && strcmp(MN(1:6),'Group_'),      % needs leading 0's
    %     fprintf('Renaming %s as ', MN);
    %     while length(MN) < 9,
    %       MN = [MN(1:6) '0' MN(7:end)];                   % pad with zeros
    %     end
    %     fprintf('%s\n', MN);
    %     FN = ['MotifLibrary' filesep MN '.mat'];
    %     Filenames(m).name = MN;
    %     save(FN,'Search','-mat');
    %   end

    fprintf('\nCurrent model is : %s\n', MN);

    Search.origmatfilename = [MN '.mat'];

    % ----- Calculate coplanar measure for each File in Search

    if ~isfield(Search.File(1),'Coplanar'),
        [L,N] = size(Search.Candidates);        % L = num instances; N = num NT
        N = N - 1;                              % number of nucleotides
        clear NewFile
        for ff = 1:length(Search.File),
            F = Search.File(ff);
            if ~isempty(F.NT),
                F.Coplanar = sparse(F.NumNT,F.NumNT);
                NewFile(ff) = F;
            end
        end
        Search.File = NewFile;

        for c = 1:length(Search.Candidates(:,1)),
            ff = Search.Candidates(c,N+1);
            i = Search.Candidates(c,1:N);

            for a = 1:length(i),
                for b = (a+1):length(i),
                    if Search.File(ff).Edge(i(a),i(b)) ~= 0,
                        NT1 = Search.File(ff).NT(i(a));
                        NT2 = Search.File(ff).NT(i(b));
                        Pair = zAnalyzePair(NT1,NT2,CL,Exemplar);
                        Search.File(ff).Coplanar(i(a),i(b)) = Pair.Coplanar;
                        Search.File(ff).Coplanar(i(b),i(a)) = Pair.Coplanar;
                    end
                end
            end
        end
        save(FN,'Search','-mat');
    end

    % ----- Calculate loose coplanar measure for each File in Search

    if ~isfield(Search.File(1),'LooseCoplanar'),
        [L,N] = size(Search.Candidates);        % L = num instances; N = num NT
        N = N - 1;                              % number of nucleotides
        clear NewFile
        for ff = 1:length(Search.File),
            F = Search.File(ff);
            if ~isempty(F.NT),
                F.LooseCoplanar = sparse(F.NumNT,F.NumNT);
                NewFile(ff) = F;
            end
        end
        Search.File = NewFile;

        for c = 1:length(Search.Candidates(:,1)),
            ff = Search.Candidates(c,N+1);
            i = Search.Candidates(c,1:N);

            for a = 1:length(i),
                for b = (a+1):length(i),
                    if Search.File(ff).Edge(i(a),i(b)) ~= 0,
                        NT1 = Search.File(ff).NT(i(a));
                        NT2 = Search.File(ff).NT(i(b));
                        Pair = zLooseCoplanar(NT1,NT2,CL,Exemplar);
                        Search.File(ff).LooseCoplanar(i(a),i(b)) = Pair.Coplanar;
                        Search.File(ff).LooseCoplanar(i(b),i(a)) = Pair.Coplanar;
                    end
                end
            end
        end
        save(FN,'Search','-mat');
    end

    % --------------------------------------- Notify

    if Verbose > 0,
        fprintf('pMakeModelsFromLibrary:  Making a JAR3D SCFG/MRF model for %s\n', MN);
    end

    % --------------------------------------- Make model and write it
    [L,N] = size(Search.Candidates);        % L = num instances; N = num NT
    N = N - 1;                              % number of nucleotides
    if L > 1 && str2num(MN(4:6)) < 900,
        % --------------------------------------- Write sequences in FASTA format
        Text = xFASTACandidates(Search.File,Search,0,MN);

        fid = fopen(['Sequences' filesep 'temp' '.fasta'],'w');
        for t = 1:length(Text),
            fprintf(fid,'%s\n',Text{t});
        end
        fclose(fid);

        sData = Alignment.loadFasta('temp.fasta');
        jF = java.lang.Boolean.FALSE;
        if strcmp(loopType,'IL'),
            seqComp = SimpleAlign.calcILEditDistances(sData,sData,jF.booleanValue());
        elseif strcmp(loopType,'HL'),
            seqComp = SimpleAlign.calcHLEditDistances(sData,sData);
        end
        if max(seqComp) > 0,
            Used = zeros(L,1);
            for i = 1:L,
                if Used(i) == 0,
                    index = 1:L;
                    index(i) = [];
                    Used(i) = 1;
                    matches = [];
                    for j = 1:length(index),
                        if seqComp(i,index(j)) == 0,
                            matches(length(matches)+1) = j;
                        end
                    end
                    index(matches) = [];
                    [Node,Truncate,Signature,RSignature] = pMakeMotifModelFromSSF(Search,Param,Prior,loopType,index);

                    Search.Signature = Signature;
                    Search.RSignature = RSignature;           %
                    save(FN,'Search','-mat');

                    MNSig = [MN '_' Signature];

                    Filenames(m).namesig = MNSig;

                    if ~strcmp(MN(4:5),'99'),                 % don't do this for helices!
                        for n = 1:length(Node),
                            if fix(abs(Node(n).Edge)) == 1,       % cWW basepair
                                %        Node(n).SubsProb = ones(4,4)/16;    % make cWW pairs noninformative
                                %        Node(n).SubsProb = [0 0 0 1; 0 0 1 0; 0 1 0 1; 1 0 1 0]/6;
                                % make cWW pairs noninformative
                            end
                        end
                    end

                    Search.Truncate = Truncate;

                    fprintf('pMakeModelsFromLibrary:  First sequence is %s\n',Text{2});

                    if ~(exist('Sequences') == 7),     % if directory doesn't yet exist
                        mkdir('Sequences');
                    end

                    if ~(exist('Models') == 7),        % if directory doesn't yet exist
                        mkdir('Models');
                    end

                    fid = fopen(['Sequences' filesep MNSig '.fasta'],'w');
                    for t = 1:length(Text),
                        fprintf(fid,'%s\n',Text{t});
                    end
                    fclose(fid);
                    Seq = strrep(Text{i*2},'*','-');
                    pWriteJavaNodeFile(Search.Query,Node,5,['L1O_' MNSig '_' Seq '.txt']);

                    % Write fasta file with omitted sequence
                    fid = fopen(['Sequences' filesep 'temp' '.fasta'],'w');
                    for t = i*2-1:i*2,
                        fprintf(fid,'%s\n',Text{t});
                    end
                    fclose(fid);

                    % Score fasta file
                    score = JAR3DMatlab.MotifParseSingle(pwd,'temp.fasta','temp.txt');
                    % quant = getQuantilesB(score,MN(4:6),loopType);
                    % Record information
                    tableIndex = tableIndex+1;
                    group(tableIndex) = {MN};
                    sequence(tableIndex) = {Text(i*2)};
                    scores(tableIndex) = score;
                    % quantiles(tableIndex) = quant;
                end

            end
        end
    end
end

% Save information
outfilename = [loopType '_L1O_diagnostic_' date '.mat'];
save(outfilename,'group','sequence','scores');