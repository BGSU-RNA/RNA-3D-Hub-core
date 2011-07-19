%----------------------------------------------------------------
function Node=GetLetters(File,Node);

m=min(min(Node.IBases));
M=max(max(Node.IBases));
Node.LeftLetter=cat(1,File.NT(Node.Left+m-1).Base);
sort(M-Node.Right(1)+Node.Right,'descend')
Node.RightLetter=cat(1,File.NT(sort(M-Node.Right(1)+Node.Right,'descend')).Base);


