
function [Seq,Scores] = makeModelScoreEmpDist(MN,loopType,sampsize,top5)

if nargin < 4,
  top5 = 0;                                  % don't analyze top 5%
end

if ~(exist(['Models' filesep 'Emp. Distributions']) == 7),
  mkdir(['Models' filesep 'Emp. Distributions']);
end

    shortName = MN(1:6);
    
    % ---------------------------------------- sample from length dist
    FN = ['Models' filesep 'Length Distributions' filesep shortName '.mat'];
    load(FN,'D','-mat');

    % ---------------------------------------- generate strand lengths

    if strcmp(loopType,'IL'),
        ld = sum(D,2);                       % marginal distribution of left
        dim=zeros(2,sampsize);
        for i = 1:sampsize,
            dim(1,i) = randsample(length(ld),1,true,ld)-1;
            rd = D(dim(1,i)+1,:);            % almost the conditional distn
            rd = rd / sum(rd);               % normalize
            dim(2,i) = randsample(length(rd),1,true,rd)-1;
        end
    elseif strcmp(loopType,'HL'),
        dim=zeros(1,sampsize);
        for i = 1:sampsize,
            dim(i) = randsample(length(D),1,true,D)-1;
        end
    end
    oFN = ['Sequences' filesep 'RS_' shortName '.fasta'];
    
    % -------------------------------------- make sequence file

    Seq = makeRandomSequencesDLFasta(oFN,dim);     % generate random sequences
    fprintf('File Written:%s\n',oFN);

    % -------------------------------------- parse and save scores

    MFN = [MN '.txt'];
    SFN = ['RS_' shortName '.fasta'];
    Scores = JAR3DMatlab.MotifParseSingle(pwd,SFN,MFN);
    distFN = ['Models' filesep 'Emp. Distributions' filesep MN(1:6) '.txt'];
    %  save(OFN,'Scores','-mat')
    fid = fopen(distFN,'w');
    for i = 1:length(Scores)
        fprintf(fid, '%f\n', Scores(i));
    end
    fclose(fid);

    % --------------------------------- delete large and cumbersome RS file

    delete(oFN);

    % --------------------------------- identify 90th percentile sequences

    y = quantile(Scores,0.90);              % top 90th percentile

fprintf('90th percentile score is %9.4f\n', y);

    i = find(Scores >= y);                  % locations of these sequences

fprintf('%d sequences score at or above the 90th percentile\n', length(i));

    MyScores = Scores(i);                   % keep just these scores
    oFN = ['Sequences' filesep 'Top10_' shortName '.fasta'];
    fid = fopen(oFN,'w');
    for k = 1:length(i),
      fprintf(fid,'> Top10 %d\n', k);
      fprintf(fid,'%s\n', Seq{i(k)});
    end
    fclose(fid);

end
