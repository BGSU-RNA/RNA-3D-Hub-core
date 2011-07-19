
% pConstructNodes
function [Node] = pConstructNodes(File,NTNumber)

[Node] = pShowStructure(File,NTNumber)
[Node] = pFineTuneNodes(File,Node);

% for n=1:length(Node)
%     [a,b]=size(Node(n).IBases);
%     if a==1
%         if abs(diff(Node(n).IBases))<5
%               Node(n).type      = 'Hairpin';                % node type
%               Node(n).subtype   = '';                   % GNRA loop
%               Node(n).nextnode  = Inf;                      % index of next node in tree
%               Node(n).P         = ones(17,1);               % state to state transitions
%               Node(n).PIns      = 1;                        % when no previous state
%               Node(n).Bm          = Node(n).IBases;
%         end
%     end
% end



OUT   = strcat(['Nodes_' File.Filename '.m']);
fidOUT      = fopen(OUT,'w+');
fprintf(fidOUT,'%% Nodes for %s', File.Filename);
fprintf(fidOUT,'%s\n',' created by pConstructNodes');

N=length(Node);
for n=1:N
    if Node(n).lpar==0
        Node(n).lpar=.001;
    end
    if Node(n).rpar==0
        Node(n).rpar=.001;
    end
      switch Node(n).type,
        case 'Initial'
            if n==1
                fprintf(fidOUT,'%s\n','n=1;  % n=1;');
            else
                fprintf(fidOUT,'%s%2d\n','n=n+1;  % n=',n);
            end
            fprintf(fidOUT,'%s%s%s\n','    Node(n).type          = ''',Node(n).type,'''');
            fprintf(fidOUT,'%s\n','    Node(n).nextnode      = n+1;');
            fprintf(fidOUT,'    Node(n).lpar          = %3g\n',Node(n).lpar(1));
            fprintf(fidOUT,'    Node(n).rpar          = %3g\n',Node(n).rpar(1));
            % tempory
            if n==1
            fprintf(fidOUT,'%s%3d%s\n','    Node(n).Bl            = [',Node(n).IBases(1),'];');
            fprintf(fidOUT,'%s%3d%s\n','    Node(n).Br            = [',Node(n).IBases(2),'];');
            end
            

        case 'Hairpin' % done
            fprintf(fidOUT,'%s%2d\n','n=n+1;  % n=',n);
            fprintf(fidOUT,'%s%s%s\n','    Node(n).type          = ''',Node(n).type,'''');
            fprintf(fidOUT,'%s\n','    Node(n).nextnode      = Inf;');
            fprintf(fidOUT,'%s\n','    Node(n).P             = ones(17,1);');
            fprintf(fidOUT,'%s\n','    Node(n).PIns          = 1;');
            fprintf(fidOUT,'%s%3d%3d%s\n','    Node(n).Bm            = [',int16(Node(n).Bm),'];');

        case 'Basepair' %done
            fprintf(fidOUT,'%s%2d\n','n=n+1;  % n=',n);
            fprintf(fidOUT,'%s%s%s\n','    Node(n).type          = ''',Node(n).type,'''');
            fprintf(fidOUT,'%s\n','    Node(n).nextnode      = n+1;');
            fprintf(fidOUT,'    Node(n).LeftLetter    = ''%s''', Node(n).LeftLetter);fprintf(fidOUT,'\n');
            fprintf(fidOUT,'    Node(n).RightLetter   = ''%s''', Node(n).RightLetter);fprintf(fidOUT,'\n');
            fprintf(fidOUT,'    Node(n).Inter         = %2d\n',Node(n).Inter);
            fprintf(fidOUT,'    Node(n).Delete        = %3g\n',Node(n).Delete);
            fprintf(fidOUT,'%s%3g%s\n','    Node(n).lpar          = [',Node(n).lpar(1),'*ones(16,1); 0];');
            fprintf(fidOUT,'%s%3g%s\n','    Node(n).rpar          = [',Node(n).rpar(1),'*ones(16,1); 0];');
            fprintf(fidOUT,'%s%3d%s\n','    Node(n).Bl            = [',int16(Node(n).Bl),'];');
            fprintf(fidOUT,'%s%3d%s\n','    Node(n).Br            = [',int16(Node(n).Br),'];');

        case 'Junction'
            fprintf(fidOUT,'%s%2d\n','n=n+1;  % n=',n);
            fprintf(fidOUT,'%s%s%s\n','    Node(n).type          = ''',Node(n).type,'''');
        case 'JunctionMotif'  
            Node(Node(n).nextnode(1)).type= 'Initial';
            Node(Node(n).nextnode(2)).type= 'Initial';
            fprintf(fidOUT,'%s%2d\n','n=n+1;  % n=',n);
            fprintf(fidOUT,'%s%s%s\n','    Node(n).type          = ''',Node(n).type,'''');
            fprintf(fidOUT,'%s%3d%4d%s\n','    Node(n).nextnode      = [',Node(n).nextnode,'];');
            fprintf(fidOUT,'%s\n','    Node(n).P             = [0.05*ones(17,1) 0.95*ones(17,1)];');
            fprintf(fidOUT,'%s\n','    Node(n).PIns          = [0.05 0.95];');
%             for t=1:length(Node(n).IBases(:,1));
%                 fprintf(fidOUT,'%s%1d%s%3d%4d%s\n','    Node(n).IBases(',t,',:)   = [',Node(n).IBases(t,:),'];');
%                 fprintf(fidOUT,'%s%d%s%2d\n','    Node(n).Inter(',t,')         = %2d\n',Node(n).Inter(t));
%             end      
%                Node(n).type      = 'JunctionMotif';               % node type
%                Node(n).nextnode  =  [a b];                   % index of next node in tree
%                Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
%                                                             % state to state transitions
%                Node(n).PIns	        = [0.05 0.95];  
            
        case 'Motif'
            fprintf(fidOUT,'%s%2d\n','n=n+1;  % n=',n);
            fprintf(fidOUT,'%s%s%s\n','    Node(n).type          = ''',Node(n).type,'''');
            fprintf(fidOUT,'%s\n','    Node(n).nextnode      = n+1;');
            fprintf(fidOUT,'%s\n','    Node(n).P             = [0.05*ones(17,1) 0.95*ones(17,1)];');
            fprintf(fidOUT,'%s\n','    Node(n).PIns          = [0.05 0.95];');
%             
%             Node(n)=BlBrLeftRight(Node(n));
%             Node(n)=GetLetters(File,Node(n));
%             Node(n).IBases=NumberIBases(Node(n).IBases);
            
            
            KL=length(Node(n).Left(1,:));
            KR=length(Node(n).Right(1,:));
            for t=1:length(Node(n).Left(:,1))
                fprintf(fidOUT,'%s%1d%s','    Node(n).Left(',t,',:)     = [');
                for k=1:KL
                    fprintf(fidOUT,'%3d',Node(n).Left(t,k));
                end
                fprintf(fidOUT,'%s\n','];');
            end
            fprintf(fidOUT,'%s\n','    Node(n).LIP           = [  1 ];'); 
            for t=1:length(Node(n).Right(:,1))
                fprintf(fidOUT,'%s%1d%s','    Node(n).Right(',t,',:)    = [');
                for k=1:KR
                    fprintf(fidOUT,'%3d',Node(n).Right(t,k));
                end
                fprintf(fidOUT,'%s\n','];');
            end
            fprintf(fidOUT,'%s\n','    Node(n).RIP           = [  1 ];'); 
            
            fprintf(fidOUT,'%s%s\n','    Node(n).LeftLetter    = ',Node(n).LeftLetter);
            fprintf(fidOUT,'%s%s\n','    Node(n).RightLetter   = ',Node(n).RightLetter);
            
            for t=1:length(Node(n).IBases(:,1));
                fprintf(fidOUT,'%s%1d%s%3d%4d%s\n','    Node(n).IBases(',t,',:)   = [',Node(n).IBases(t,:),'];');
                fprintf(fidOUT,'%s%d%s%3g\n','    Node(n).Inter(',t,')      = ',Node(n).Inter(t));
            end
            
            if diff(Node(n).Bl)==0 
                fprintf(fidOUT,'%s%3d%s\n','    Node(n).Bl            = [',Node(n).Bl(1),'];');
            else
                fprintf(fidOUT,'%s%3d%s%3d%s\n','    Node(n).Bl            = [',Node(n).Bl(1),',',Node(n).Bl(2),'];');
            end
            if diff(Node(n).Br)==0 
                fprintf(fidOUT,'%s%3d%s\n','    Node(n).Br            = [',Node(n).Br(1),'];');
            else
                fprintf(fidOUT,'%s%3d%s%3d%s\n','    Node(n).Bl            = [',Node(n).Br(1),',',Node(n).Br(2),'];');
            end
            
            
%               Node(n).Left(1,:) = [1 2 3];                  % nucleotides to use on left
%               Node(n).Left(2,:) = [1 3 4];                  % nucleotides to use on left
%               Node(n).Left(3,:) = [1 2 4];                  % nucleotides to use on left
%               Node(n).LIP       = [0.8 0.1 0.1];            % probs for insertion possibs
%               Node(n).Right(1,:)= [4 3 2 1];                % nucleotides to use on right
%               Node(n).RIP       = 1;                        % probs for insertion possibs

        case 'Alternative'
            fprintf(fidOUT,'%s%2d\n','n=n+1;  % n=',n);
      end
      fprintf(fidOUT,'\n');
end


fclose(fidOUT);



% %----------------------------------------------------------------
% function Node=BlBrLeftRight(Node)
% 
% L=length(Node.IBases(:,1))*2;
% Bases=reshape(Node.IBases,1,L);
% D     = diff(Bases);
% [a,b] = max(D);
% 
% % Bl and Br
% Node.Bl=[min(Bases(1:b)),max(Bases(1:b))];
% Node.Br=[min(Bases(b+1:L)),max(Bases(b+1:L))];
% 
% % Left and Right
% Node.Left(1)=1;
% B1=unique(Bases(1:b));
% B2=unique(Bases(b+1:L));
% D1=diff(B1);
% D2=diff(B2);
% for i=2:length(B1)
%     Node.Left(i)=Node.Left(i-1)+D1(i-1);
% end
% RightTemp(1)=1;
% for i=2:length(B2)
%     RightTemp(i)=RightTemp(i-1)+D2(i-1);
% end
% Node.Right=sort(RightTemp,'descend');
% 
% 
% %----------------------------------------------------------------
% function Node=GetLetters(File,Node);
% 
% m=min(min(Node.IBases));
% M=max(max(Node.IBases));
% Node.LeftLetter=cat(1,File.NT(Node.Left+m-1).Base);
% Node.RightLetter=cat(1,File.NT(sort(M-Node.Right(1)+Node.Right,'descend')).Base);
% 
% 
% %----------------------------------------------------------------
% function IBases=NumberIBases(IBases)
% 
% L = length(IBases(:,1));
% [a,b] = unique(IBases);
% for k=1:length(a)
%    c= find(IBases==a(k));
%    for t=1:length(c)
%         IBases(c)=k;
%    end
% end
% 
% 
