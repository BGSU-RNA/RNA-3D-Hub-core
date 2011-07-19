
% Example:
% BetterEmpDist('IL_225_05_cWW-cWW-cSH','IL',100)

function [Seq,Scores] = BetterEmpDist(MN,loopType,sampsize,prec,top10)

if nargin < 4,
    prec = 4;
end

if nargin < 5,
  top10 = 0;                             % don't analyze top 10% of sequences
end

if ~(exist(['Models' filesep 'Emp. Distributions']) == 7),
    mkdir(['Models' filesep 'Emp. Distributions']);
end

shortName = MN(1:6);

% ---------------------------------------- load length dist
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

FN
D

    dim=zeros(1,sampsize);
    for i = 1:sampsize,
        dim(i) = randsample(length(D),1,true,D)-1;
    end
end
oFN = ['Sequences' filesep 'RS_' shortName '.fasta'];

% -------------------------------------- make sequence file

Seq = makeRandomSequencesDLFasta(oFN,dim);     % generate random sequences
fprintf('File written:%s\n',oFN);

% -------------------------------------- parse

MFN = [MN '.txt'];
SFN = ['RS_' shortName '.fasta'];
Scores = JAR3DMatlab.MotifParseSingle(pwd,SFN,MFN);

% -------------------------------------- calculate distribution
Values=unique(Scores);
N=histc(Scores,Values);
step=1/length(Scores);
dist = zeros(length(Values),1);
dist(1) = N(1)*step;
for i = 2:length(Values)
    dist(i) = dist(i-1) + N(i) * step;
end

% -------------------------------------- save distribution
distFN = ['Models' filesep 'Emp. Distributions' filesep MN(1:6) '.txt'];
fid = fopen(distFN,'w');
dist=round(dist*(10^prec))/(10^prec);
prev = -10000;
for i = 1:length(dist)
    cur = dist(i);
    if cur > prev,
        fprintf(fid,'%f %.4f\n',Values(i),dist(i));
    end
    prev = cur;
end
fclose(fid);

% --------------------------------------- delete large and cumbersome RS file

delete(oFN);

% --------------------------------- identify 90th percentile sequences

if top10 > 0,
  y = quantile(Scores,0.90);              % top 90th percentile
  fprintf('90th percentile score is %9.4f\n', y);

  i = find(Scores >= y);                  % locations of these sequences

  fprintf('%d sequences score at or above the 90th percentile\n', length(i));

  oFN = ['Sequences' filesep 'Top10_' shortName '.fasta'];
  fid = fopen(oFN,'w');
  for k = 1:length(i),
    fprintf(fid,'> Top10 %d\n', k);
    fprintf(fid,'%s\n', Seq{i(k)});
  end
  fclose(fid);
end
