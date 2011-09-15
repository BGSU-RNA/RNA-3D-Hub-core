function NormC = pClusterNorm(InterIndicesF,SubsProbF,LeftIndex,RightIndex)
    psum = 0;
    [numInterF,dum] = size(InterIndicesF);
    InterIndices = InterIndicesF(1,:);
    Inters = 1;
    added = 1;
    while added == 1,
        added = 0;
        for i = 2:numInterF,
            if ~isempty(intersect(InterIndices(:),InterIndicesF(i,:))) && isempty(intersect(i,Inters)),
                [numInter,dum] = size(InterIndices);
                InterIndices(numInter+1,:) = InterIndicesF(i,:);
                Inters(numInter+1) = i;
                added = 1;
            end
        end
    end
    SubsProb = SubsProbF(:,:,Inters);
    IBases = sort(unique(InterIndices(:)));
    Left = [];  Right = [];
    for i = 1:length(IBases),
        if ismember(IBases(i),LeftIndex);
            Left(length(Left)+1) = find(LeftIndex == IBases(i));
        elseif ismember(IBases(i),RightIndex);
            Right(length(Right)+1) = find(RightIndex == IBases(i));
        end
    end
    numBases = length(Left) + length(Right);
    [numInter,dum] = size(InterIndices);
    Cols(1:length(Left)) = LeftIndex(Left);
    Cols(length(Left)+1:length(Left)+length(Right)) = RightIndex(Right);
    % System.out.println("ClusterNode.Normalize "+numIndices);
    if numInter == 1,
        psum = 1;
    else
        code = ones(1,numBases);
        for i = 1:4^numBases,
            [numInter,dum] = size(InterIndices);
            % score codes[] according to the various interactions
            prob = 0;
            for j = 1:numInter,
                i1 = find(Cols == InterIndices(j,1));
                i2 = find(Cols == InterIndices(j,2));
                prob = prob+log(SubsProb(code(i1),code(i2),j));
            end
            increm = 0;
            j = 1;
            if i < 4^numBases,    %if we're not done yet, increment BP code
                while increm == 0,
                    if code(j) == 4,
                        code(j) = 1;
                        j = j + 1;
                    else
                        code(j) = code(j) +1;
                        increm = 1;
                    end
                end
            end
            psum = psum+exp(prob);
        end
    end
    if numInterF-numInter <= 1,
        NormC = psum;
    else
        remInt = setDiff(1:numInterF,Inters);
        InterIndicesR = InterIndicesF(remInt,:);
        SubsProbR = SubsProbF(:,:,remInt);
        NormC = psum*pClusterNorm(InterIndicesR,SubsProbR,LeftIndex,RightIndex);
    end
end