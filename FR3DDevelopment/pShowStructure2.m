% pShowStructure2(File,NTNumber) displays the secondary structure of a
% molecule as reflected by the interaction matrix, starting at NTNumber and
% ending with whather base NTNumber interacts with

function [Node] = pShowStructure2(File,NTNumber,LastNTNumber,Interact,Node,n)

% if File is a text string (filename), load the file and display
if nargin<6
   n=0;
end
if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

% if NTNumber is a cell array of numbers, look up the indices

if strcmp(class(NTNumber),'char'),
  NTNumber = {NTNumber};
end

if strcmp(class(NTNumber),'cell'),
  Indices = zIndexLookup(File,NTNumber);
else
  Indices = NTNumber;
end

N = length(File.NT);                           % number of nucleotides
G = (File.Inter < 15) .* (File.Inter ~= 0);    % consider basepairing only

if nargin == 2,
  LastNTNumber = N;
elseif strcmp(class(LastNTNumber),'cell'),
  LastNTNumber = zIndexLookup(File,LastNTNumber);
end

if nargin < 4,
 for a = 1:N,                              % loop through nucleotides
  k = find(G(a,:));                        % find indices of interacting bases
  [y,L] = sort(abs(File.Inter(a,k)));      % sort by interaction category
  Interact{a}.Categ = File.Inter(a,k(L));  % store categories
  Interact{a}.Index = k(L);                % store indices of interacting bases
 end
end

a = Indices(1);                            % first index
b = Inf;                                   % index of base it interacts with
c = 2;                                     % category of interaction
t = 0;
A  = a;
B  = LastNTNumber;
AA = a;
BB = a;
% Initial node creation ------
n=n+1;
Node(n).type          = 'Initial';                % node type
Node(n).nextnode      = n+1;                      % index of next node in tree
Node(n).lpar          = 0;                        % left insertion parameter
Node(n).rpar          = 0; 
Node(n).Bl            = 0;
Node(n).Br            = 0;
I=1;
% -----------------------------
BP=0;


while ((a < b) | (abs(c) > 1)) & (a <= N), % while not the end of the loop,

  if length(Interact{a}.Categ) > 0,        % if a interacts with something,
    b = Interact{a}.Index(1);              % get index 
    c = Interact{a}.Categ(1);              % and category of interaction
    if (a < b),                            % if b comes after a

      % ---------------------------------- Check for junction

      if (sum(sum(File.Inter((b+1):(BB-1),(b+1):(BB-1)) == 1)) > 0) & ...
         (sum(sum(File.Inter((a+1):(b-1),(a+1):(b-1))   == 1)) > 0),
           I=0;
           n=n+1;
           Node(n).type      = 'JunctionMotif';               % node type
           Node(n).nextnode(1)  =  n+1;                   % index of next node in tree
           Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
                                                        % state to state transitions
           Node(n).PIns	        = [0.05 0.95];                                % when no previous state
           Node=pShowStructure2(File,a,b,Interact,Node,n);
           Node(n).nextnode(2)= length(Node)+1;
           Node=pShowStructure2(File,b+1,BB,Interact,Node,length(Node));
        return
      end

        % ----------------------------- check for pasepairs
      if ((sum(sum(File.Inter((a+1):(b-1),(a+1):(b-1)) == 1)) == 0) & ...
          (sum(sum(File.Inter((a+1):B,    (a+1):B)     == 1))  > 0)) | ...
         (b > B),
        d = 1;
      else
        % ------------------------- Show interaction between a and b
        d = 2;
        A = a;                               % a of last pair
        B = b;                               % b of last pair
        if c == 1,
          AA = a;                            % a of last class 1 pair
          BB = b;                            % b of last class 1 pair
        end
      end
    if ((length(Interact{a}.Categ) == 1) & (length(Interact{b}.Categ) == 1))  
          t=0;
          I=0;
          n=n+1;  
          BP=1;
          Node(n).type        = 'Basepair';                % node type
          Node(n).nextnode    = n+1;                      % index of next node in tree
          Node(n).LeftLetter  = File.NT(a).Base;
          Node(n).RightLetter = File.NT(b).Base;
          Node(n).Inter       = c;
          Node(n).Delete      = 0.1;
          Node(n).lpar        = [0.001*ones(16,1); 0];     % left insertion parameters
          Node(n).rpar        = [0.001*ones(16,1); 0];     % right insertion parameters
          Node(n).Bl          = [str2num(File.NT(a).Number)];
          Node(n).Br          = [str2num(File.NT(b).Number)];  
          Node(n).Inter       = c;
          Node(n).IBases      = [str2num(File.NT(a).Number) str2num(File.NT(b).Number)];   
    elseif ((length(Interact{a}.Categ) >1 ) & (length(Interact{b}.Categ) >= 1)) | ...
            ((length(Interact{a}.Categ) >=1 ) & (length(Interact{b}.Categ) > 1))
          n=n+1;
          t=0;
          Node(n).type = 'Motif';
          BP=0;
          Node(n).lpar        = [0.001*ones(16,1); 0];     % left insertion parameters
          Node(n).rpar        = [0.001*ones(16,1); 0];     % right insertion parameters 
        for k = 1:length(Interact{a}.Categ),
          bb = Interact{a}.Index(k);
          cc = Interact{a}.Categ(k);
          t=t+1;
          Node(n).IBases(t,:) = [str2num(File.NT(a).Number) str2num(File.NT(bb).Number)];
          Node(n).Inter(t)    = cc;
        end
          Node(n).type = 'Motif';
          BP=0;
          Node(n).lpar        = [0.001*ones(16,1); 0];     % left insertion parameters
          Node(n).rpar        = [0.001*ones(16,1); 0];     % right insertion parameters
      for k=1:length(Interact{b}.Categ),
        bb = Interact{b}.Index(k);
        cc = Interact{b}.Categ(k);
        if bb ~= a,
          t=t+1;
          if b<bb
             Node(n).IBases(t,:) = [str2num(File.NT(b).Number) str2num(File.NT(bb).Number)];
          else
             Node(n).IBases(t,:) = [str2num(File.NT(bb).Number) str2num(File.NT(b).Number)]; 
          end
          Node(n).Inter(t)    = cc;           
        end
      end
      
    else
        
        
    end
    end
  end

a = a + 1;
end

