

disp('Make sure the Matlab current folder is My Dropbox\BGSURNA\Motifs');
sampsize = 100000;
loopType = 'HL';                            % hairpins only
loopType = 'IL';                            % internal loops only

Filenames = dir(['Models' filesep 'Length Distributions' filesep]);

for m = 1:length(Filenames),
    if (length(Filenames(m).name) > 2),
      if strcmp(Filenames(m).name(1:2),loopType), 
          keep(m) = 1;
          Filenames(m).name = strrep(Filenames(m).name,'.mat','');
      end
    end 
end
Filenames = Filenames(find(keep));

for m = 1:length(Filenames),
    MN = Filenames(m).name;
    FN = ['Models' filesep 'Length Distributions' filesep MN '.mat'];
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
            dim(i,1) = randsample(length(D),1,true,D)-1;
        end
    end
    oFN = ['Sequences' filesep 'RS_' MN '.fasta'];

    makeRandomSequencesDLFasta(oFN,dim);     % generate random sequences
    fprintf('File Written:%s\n',oFN);
    
end

pEstModelScoreDistr
