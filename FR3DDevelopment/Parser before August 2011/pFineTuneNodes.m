%pFineTuneNodes

Node2 = Node;

%function [Node] = pFineTuneNodes(File,Node2);

N=length(Node2);
T=0;
count=0;
Junction=0;
for n=1:N
    T=T+1;
    Node2(n).Left=[];
    Node2(n).Right=[];
    if Node2(n).lpar==0
        Node2(n).lpar=.001;
    end
    if Node2(n).rpar==0
        Node2(n).rpar=.001;
    end
      switch Node2(n).type,
        case 'Initial'
                    Node(T).type        = Node2(n).type;
                    Node(T).nextnode    = Node2(n).nextnode;
                    Node(T).lpar        = Node2(n).lpar;
                    Node(T).rpar        = Node2(n).rpar;
                    Node(T).Bl          = Node2(n).Bl;
                    Node(T).Br          = Node2(n).Br;
                    Node(T).IBases      = Node2(n).IBases;
                    Node(T).Edge       = Node2(n).Edge;
                    Node(T).LeftLetter  = Node2(n).LeftLetter;
                    Node(T).RightLetter = Node2(n).RightLetter;
                    Node(T).Delete      = Node2(n).Delete;
                    Node(T).P           = Node2(n).P;
                    Node(T).PIns        = Node2(n).PIns;
                    Node(T).Left        = Node2(n).Left;
                    Node(T).Right       = Node2(n).Right;

        case 'Hairpin'
                    Node(T).type        = Node2(n).type;
                    Node(T).nextnode    = Node2(n).nextnode;
                    Node(T).lpar        = Node2(n).lpar;
                    Node(T).rpar        = Node2(n).rpar;
                    Node(T).Bl          = Node2(n).Bl;
                    Node(T).Br          = Node2(n).Br;
                    Node(T).IBases      = Node2(n).IBases;
                    Node(T).Edge       = Node2(n).Edge;
                    Node(T).LeftLetter  = Node2(n).LeftLetter;
                    Node(T).RightLetter = Node2(n).RightLetter;
                    Node(T).Delete      = Node2(n).Delete;
                    Node(T).P           = Node2(n).P;
                    Node(T).PIns        = Node2(n).PIns;
                    Node(T).Left        = Node2(n).Left;
                    Node(T).Right       = Node2(n).Right;
                    Node(T).NL          = Node2(n).NL;
                    Node(T).NR          = Node2(n).NR;
                    Node(T).IBasesOld   = Node2(n).IBasesOld;
        case 'Basepair'
                    Node(T).type        = Node2(n).type;
                    Node(T).nextnode    = Node2(n).nextnode;
                    Node(T).lpar        = Node2(n).lpar;
                    Node(T).rpar        = Node2(n).rpar;
                    Node(T).Bl          = Node2(n).Bl;
                    Node(T).Br          = Node2(n).Br;
                    Node(T).IBases      = Node2(n).IBases;
                    Node(T).Edge       = Node2(n).Edge;
                    Node(T).LeftLetter  = Node2(n).LeftLetter;
                    Node(T).RightLetter = Node2(n).RightLetter;
                    Node(T).Delete      = Node2(n).Delete;
                    Node(T).P           = Node2(n).P;
                    Node(T).PIns        = Node2(n).PIns;
                    Node(T).Left        = Node2(n).Left;
                    Node(T).Right       = Node2(n).Right;

        case 'Junction'
                    Node(T).type        = Node2(n).type;
                    Node(T).nextnode    = Node2(n).nextnode;
                    Node(T).lpar        = Node2(n).lpar;
                    Node(T).rpar        = Node2(n).rpar;
                    Node(T).Bl          = Node2(n).Bl;
                    Node(T).Br          = Node2(n).Br;
                    Node(T).IBases      = Node2(n).IBases;
                    Node(T).Edge       = Node2(n).Edge;
                    Node(T).LeftLetter  = Node2(n).LeftLetter;
                    Node(T).RightLetter = Node2(n).RightLetter;
                    Node(T).Delete      = Node2(n).Delete;
                    Node(T).P           = Node2(n).P;
                    Node(T).PIns        = Node2(n).PIns;
                    Node(T).Left        = Node2(n).Left;
                    Node(T).Right       = Node2(n).Right;
                    Node(T).NL          = Node2(n).NL;
                    Node(T).NR          = Node2(n).NR;
                    Node(T).IBasesOld   = Node2(n).IBasesOld;
                    
        case 'JunctionMotif'  
                    Junction=T;
                    Node(T).type        = Node2(n).type;
                    Node(T).nextnode    = Node2(n).nextnode;
                    Node(T).lpar        = Node2(n).lpar;
                    Node(T).rpar        = Node2(n).rpar;
                    Node(T).Bl          = Node2(n).Bl;
                    Node(T).Br          = Node2(n).Br;
                    Node(T).IBases      = Node2(n).IBases;
                    Node(T).Edge       = Node2(n).Edge;
                    Node(T).LeftLetter  = Node2(n).LeftLetter;
                    Node(T).RightLetter = Node2(n).RightLetter;
                    Node(T).Delete      = Node2(n).Delete;
                    Node(T).P           = Node2(n).P;
                    Node(T).PIns        = Node2(n).PIns;
                    Node(T).Left        = Node2(n).Left;
                    Node(T).Right       = Node2(n).Right;
%                     Node(T).NL          = Node2(n).NL;
%                     Node(T).NR          = Node2(n).NR;
%                     Node(T).IBasesOld   = Node2(n).IBasesOld;

        case 'Motif'
            
            Node2(n).NL = 0;
            Node2(n).NR = 0;
            Node2(n)=BlBrLeftRight(Node2(n));

            k=length(Node2(n).IBases(:,1));
            K=ones(1,k); 
            for i=1:k-1
                for j=i+1:k
                    if K(i) & length(setdiff(Node2(n).IBases(i,:),Node2(n).IBases(j,:)))==0
                       K(j) = 0; 
                    end
                end
            end
            K=find(K); 
            IBases=[];
            Edge=[];
            for i=1:length(K)
                IBases(i,:) = sort(Node2(n).IBases(K(i),:));
                Edge(i)   = Node2(n).Edge(K(i)); 
            end
            
            Node2(n).IBases=IBases;
            Node2(n).Edge=Edge;
            [a,b] = sort(Node2(n).IBases(:,1));
            Node2(n).IBases = Node2(n).IBases(b,:);
            Node2(n).Edge  = Node2(n).Edge(b);
            Node2(n).Br     = Node2(n).Br;
            Node2(n).IBasesOld=Node2(n).IBases;
            Node2(n)=GetLetters(File,Node2(n));
            Node2(n)=NumberIBases(Node2(n));
            TempNode=CheckSplitNode(Node2(n));
n     
            P=length(TempNode);
            if P>1 % change back to 1
                for p=1:P
                   count=count+1;
                   TNextNode(count)=T;
                   TempNode(p).Bl     = [];
                   TempNode(p).Br     = [];
                   TempNode(p).Left   = [];
                   TempNode(p).NL     = [];
                   TempNode(p).NR     = [];
                   TempNode(p).Right  = [];
                   TempNode(p).IBases = TempNode(p).IBasesOld;
                   TempNode(p) = BlBrLeftRight(TempNode(p));
                   TempNode(p) = NumberIBases(TempNode(p));
                   TempNode(p).LeftLetter  = '';
                   TempNode(p).RightLetter = '';
TempNode(p).IBases
TempNode(p)
                   TempNode(p)=GetLetters(File,TempNode(p));
                    Node(T).type        = Node2(n).type;
                    Node(T).nextnode    = Node2(n).nextnode;
                    Node(T).lpar        = Node2(n).lpar;
                    Node(T).rpar        = Node2(n).rpar;
                    Node(T).Bl          = TempNode(p).Bl;
                    Node(T).Br          = TempNode(p).Br;
                    Node(T).IBases      = TempNode(p).IBases;
                     Node(T).Edge       = TempNode(p).Edge;
                    Node(T).LeftLetter  = TempNode(p).LeftLetter;
                    Node(T).RightLetter = TempNode(p).RightLetter;
                    Node(T).Delete      = Node2(n).Delete;
                    Node(T).P           = Node2(n).P;
                    Node(T).PIns        = Node2(n).PIns;
                    Node(T).Left        = TempNode(p).Left;
                    Node(T).Right       = TempNode(p).Right;
                    Node(T).NL          = TempNode(p).NL;
                    Node(T).NR          = TempNode(p).NR;
                    Node(T).IBasesOld   = TempNode(p).IBasesOld;
                    T=T+1;
                end
                T=T-1;   
            else
                    Node(T).type        = Node2(n).type;
                    Node(T).nextnode    = Node2(n).nextnode;
                    Node(T).lpar        = Node2(n).lpar;
                    Node(T).rpar        = Node2(n).rpar;
                    Node(T).Bl          = Node2(n).Bl;
                    Node(T).Br          = Node2(n).Br;
                    Node(T).IBases      = Node2(n).IBases;
                    Node(T).Edge       = Node2(n).Edge;
                    Node(T).LeftLetter  = Node2(n).LeftLetter;
                    Node(T).RightLetter = Node2(n).RightLetter;
                    Node(T).Delete      = Node2(n).Delete;
                    Node(T).P           = Node2(n).P;
                    Node(T).PIns        = Node2(n).PIns;
                    Node(T).Left        = Node2(n).Left;
                    Node(T).Right       = Node2(n).Right;
                    Node(T).NL          = Node2(n).NL;
                    Node(T).NR          = Node2(n).NR;
                    Node(T).IBasesOld   = Node2(n).IBasesOld;
            end
            
%             if strcmp(Node2(n+1).type,'Motif')
%                 TempNode=CheckNodeOverlap(Node2,n);
%             end

        case 'Alternative'
                    Node(T).type        = Node2(n).type;
                    Node(T).nextnode    = Node2(n).nextnode;
                    Node(T).lpar        = Node2(n).lpar;
                    Node(T).rpar        = Node2(n).rpar;
                    Node(T).Bl          = Node2(n).Bl;
                    Node(T).Br          = Node2(n).Br;
                    Node(T).IBases      = Node2(n).IBases;
                    Node(T).Edge       = Node2(n).Edge;
                    Node(T).LeftLetter  = Node2(n).LeftLetter;
                    Node(T).RightLetter = Node2(n).RightLetter;
                    Node(T).Delete      = Node2(n).Delete;
                    Node(T).P           = Node2(n).P;
                    Node(T).PIns        = Node2(n).PIns;
                    Node(T).Left        = Node2(n).Left;
                    Node(T).Right       = Node2(n).Right;
                    Node(T).NL          = Node2(n).NL;
                    Node(T).NR          = Node2(n).NR;
                    Node(T).IBasesOld   = Node2(n).IBasesOld;
      end
end
if Junction>0                              % if there is a junction motif
    a=find(TNextNode<Node(Junction).nextnode(1));
    b=find(TNextNode<Node(Junction).nextnode(2));
    Node(Junction).nextnode(1)=Node(Junction).nextnode(1)+length(a);
    Node(Junction).nextnode(2)=Node(Junction).nextnode(2)+length(b);
    Node(Node(Junction).nextnode(2)-1).type='Hairpin';
    Node(Node(Junction).nextnode(2)-1).Bm=[min(Node(Node(Junction).nextnode(2)-1).Bl),...
        max(Node(Node(Junction).nextnode(2)-1).Br)];
end

Node(T).type='Hairpin';
Node(T).Bm=[min(Node(T).Bl),max(Node(T).Br)];

