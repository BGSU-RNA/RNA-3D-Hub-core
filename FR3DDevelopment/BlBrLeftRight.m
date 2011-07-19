%------------------------------------------
function Node=BlBrLeftRight(Node)

L=length(Node.IBases(:,1))*2;
Bases=sort(reshape(Node.IBases,1,L));
D     = diff(Bases);
[a,b] = max(D);

% Bl and Br
Node.Bl=[min(Bases(1:b)),max(Bases(1:b))];
Node.Br=[min(Bases(b+1:L)),max(Bases(b+1:L))];

% Left and Right
Node.Left(1)=1;
B1=unique(Bases(1:b));
B2=unique(Bases(b+1:L));
Node.NL=length(B1);
Node.NR=length(B2);
D1=diff(B1);
D2=diff(B2);
for i=2:length(B1)
    Node.Left(i)=Node.Left(i-1)+D1(i-1);
end
RightTemp(1)=1;
for i=2:length(B2)
    RightTemp(i)=RightTemp(i-1)+D2(i-1);
end
Node.Right=sort(RightTemp,'descend');


