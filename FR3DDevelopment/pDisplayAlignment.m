function pDisplayAlignment(Alig,probsM,jj,graph,fileprint,filename)
  if nargin < 3
    jj = -1;
  end
  if nargin < 4
    graph = 0;
  end
  if nargin < 5 
    fileprint = 0;
  end
  if fileprint
      if nargin < 6
          error('Please specify a filename');
      end
  end
  diffMat = diff(probsM(:,:),1,2);
  sortrange = length(diffMat(1,:));
  %sortrange = 20;
  D = zMutualDistance(diffMat(:,1:sortrange),200);
  p = zOrderbySimilarity(D);
  if fileprint
    fid = fopen(filename,'w+');
  else
    fid = 1;
  end
  [NumSeq,dummy]=size(probsM);
  fprintf(fid,'JAR3D %d\n',jj);
  fprintf(fid,'%s\n', Alig.get(0));
  for i = 1:NumSeq,
      j=p(i);
      if 3*j <= size(Alig),
        fprintf(fid,'%s %30s %s\n', Alig.get(3*j-2), Alig.get(3*j-1), Alig.get(3*j));  
      end
  end
  if fileprint
      fclose(fid);
  end
  if graph
    pcolor(real(diffMat(p,:)));
    shading flat
    caxis([0 10])
    figure
    zGraphDistanceMatrix(D(p,p));
  end
end