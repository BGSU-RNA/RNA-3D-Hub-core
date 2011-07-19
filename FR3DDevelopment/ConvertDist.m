
ModelFNs = dir(['Models' filesep]);
keep = zeros(1,length(ModelFNs));
for m = 1:length(ModelFNs),
    if (length(ModelFNs(m).name) > 2),
        if (strcmp(ModelFNs(m).name(1:2),loopType) && ~strcmp(ModelFNs(m).name(4),'M') && ~strcmp(ModelFNs(m).name(4),'.')),
            keep(m) = 1;
            ModelFNs(m).name = strrep(ModelFNs(m).name,'.txt','');
        end
    end
end
ModelFNs = ModelFNs(find(keep));

for m = 1:length(ModelFNs),
    MN = ModelFNs(m).name;
    distFN = ['Models' filesep 'Emp. Distributions' filesep MN(1:6) '.txt'];
    Scores = load(distFN);
    prec = 4;
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
end
