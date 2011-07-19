%------------------------------------------
function TN=CheckSplitNode(Node);

L=length(Node.IBases(:,1));
M=Node.NL+Node.NR;
for k=1:M
    [a,b]=find(Node.IBases==k);
    Group(1,k)=2*L;
    Group(2,k)=0;
    for i=1:length(a)
        Group(1,k)=min(max(Node.IBases(a(i),mod(b(i),2)+1),k),Group(1,k));
        Group(2,k)=max(max(Node.IBases(a(i),mod(b(i),2)+1),k),Group(2,k));
    end
end

G=zeros(1,M);
g=1;
k=1;

while min(G)==0
    a=find(((G==0).*Group(2,:))>=Group(1,k));
    count=0;
    while count<2
        count=count+1;
        G(a)=g;
        [v,b]=min(Group(1,a));
        c=a(b);
        a=find(((G==0).*Group(2,:))>=Group(1,c));
        if length(a)==0
            count=count+1;
        end
    end
    g=g+1;
    [v,k]=min(G);
end
if max(G)==1
    TN=0;
else % need to split node
    TN=[];
    N=max(G);
    for k=1:N
        a=find(G==k);
        count=0;
        for i=1:L
            if length(intersect(a,Node.IBases(i,:)))>0
                count=count+1;
                TN(k).IBases(count,:)=Node.IBases(i,:);
                TN(k).IBasesOld(count,:)=Node.IBasesOld(i,:);
                TN(k).Edge(count)=Node.Edge(i);
            end  
        end            
    end
end

% Node.IBasesOld
% disp('--------------------')

% pause
