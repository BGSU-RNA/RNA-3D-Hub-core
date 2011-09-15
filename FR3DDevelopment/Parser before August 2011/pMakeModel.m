% pMakeModel(File,NTNumber) makes a secondary structure node model based on
% the Edge interaction matrix in File, starting at NTNumber and ending at
% LastNTNumber

function [Node] = pMakeModel(File,NTNumber,LastNTNumber,Interact,Node,n)

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

N = length(File.NT);                       % number of nucleotides
E = abs(fix(File.Edge));                   % don't distinguish subcategories
G = (E < 15) .* (E ~= 0);                  % consider basepairing only

if nargin == 2,
  LastNTNumber = N;
elseif strcmp(class(LastNTNumber),'cell'),
  LastNTNumber = zIndexLookup(File,LastNTNumber);
end

if nargin < 4,
 for a = 1:N,                              % loop through nucleotides
  k = find(G(a,:));                        % find indices of interacting bases
  [y,L] = sort(E(a,k));                    % sort by edge interaction category
  Interact{a}.Categ = File.Edge(a,k(L));   % store categories
  Interact{a}.Index = k(L);                % store indices of interacting bases
 end
end

a = Indices(1);                            % first index
b = Inf;                                   % index of base it interacts with
c = 2;                                     % category of interaction
t = 0;

A  = a;                                    % previous interacting base on left
B  = LastNTNumber;                         % previous interacting base on right
AA = a;                                    % previous cWW base on left
BB = a;                                    % previous cWW base on right

% Initial node creation ------
n = n+1;
Node(n).type      = 'Initial';                % node type
Node(n).nextnode  = n+1;                      % index of next node in tree
Node(n).lpar      = 0;                        % left insertion parameter
Node(n).rpar      = 0; 
Node(n).Bl        = 0;                        % actual left base
Node(n).Br        = 0;
I                 = 1;

% -----------------------------

BP=0;
while ((a < b) | (abs(c) > 1)) & (a <= N), % while not the end of the loop,

  if length(Interact{a}.Categ) > 0,        % if a interacts with something,
    b = Interact{a}.Index(1);              % get index of that something
    c = Interact{a}.Categ(1);              % and category of interaction

    if (a < b),                            % if b comes after a

      % ---------------------------------- Check for junction

      if (sum(sum(E((b+1):(BB-1),(b+1):(BB-1)) == 1)) > 0) & ...
         (sum(sum(E((a+1):(b-1),(a+1):(b-1))   == 1)) > 0),
           n = n+1;
           Node(n).type        = 'JunctionMotif';  % node type
           Node(n).nextnode(1) =  n+1;          % index of next node in tree
           Node = pMakeModel(File,a,b,Interact,Node,n);
           Node(n).nextnode(2)  = length(Node)+1;
           Node = pMakeModel(File,b+1,BB,Interact,Node,length(Node));

           Node(n).P         = [0.05*ones(17,1) 0.95*ones(17,1)];
                                                 % state to state transitions
           Node(n).PIns	        = [0.05 0.95];   % when no previous state

        return
      end

      % ---------------------------------- Check for insertions on right

      if (b < B-1),                         % gap from last interacting base
        if sum(sum(E((b+1):(B-1),(b+1):(B-1)) == 1)) == 0, % no helix there
          for e = (B-1):-1:(b+1),           % loop through insertions
            if strcmp(Node(n).type,'Initial')
              Node(n).rpar=Node(n).rpar+1;
            else
              Node(n).rpar=Node(n).rpar+[ones(16,1); 0];
            end
            if length(Interact{e}.Categ) > 0,
             for k=1:length(Interact{e}.Categ),
              bb = Interact{e}.Index(k);
              cc = Interact{e}.Categ(k);
                t=t+1;
                Node(n).type = 'Motif';
                Node(n).Edge(t)     = cc;
                Node(n).IBases(t,:) = [str2num(File.NT(e).Number) str2num(File.NT(bb).Number)];
             end
            end
          end
        end
      end

      % ---------------------------------- Check if b is out of sequence

      if ((sum(sum(E((a+1):(b-1),(a+1):(b-1)) == 1)) == 0) & ...
          (sum(sum(E((a+1):B,    (a+1):B)     == 1))  > 0)) | ...
         (b > B),
        d = 1;
      else
        % ---------------------------------- Show interaction between a and b
    if ((length(Interact{a}.Categ) == 1) & (length(Interact{b}.Categ) == 1))  
          t=0;
          n=n+1;  
          BP=1;
          Node(n).type        = 'Basepair';                % node type
          Node(n).nextnode    = n+1;                      % index of next node in tree
          Node(n).LeftLetter  = File.NT(a).Base;
          Node(n).RightLetter = File.NT(b).Base;
          Node(n).Edge       = c;
          Node(n).Delete      = 0.1;
          Node(n).lpar        = [.001*ones(16,1); 0];     % left insertion parameters
          Node(n).rpar        = [.001*ones(16,1); 0];     % right insertion parameters
          Node(n).Bl          = [str2num(File.NT(a).Number)]; %nucleotide number
          Node(n).Br          = [str2num(File.NT(b).Number)]; % nucl number 
          Node(n).Edge       = c;
          Node(n).IBases      = [str2num(File.NT(a).Number) str2num(File.NT(b).Number)];   
    else
        t=t+1;
        Node(n).IBases(t,:)   = [str2num(File.NT(a).Number) str2num(File.NT(b).Number)]; 
        Node(n).Edge(t)      = c;
    end
        d = 2;
        A = a;                               % a of last pair
        B = b;                               % b of last pair
        if c == 1,
          AA = a;                            % a of last class 1 pair
          BB = b;                            % b of last class 1 pair
        end
      end
    elseif (c ~= 1) & (a < B),
        if strcmp(Node(n).type,'Initial')
            Node(n).rpar=Node(n).rpar+1;
        else
            Node(n).rpar=Node(n).rpar+[ones(16,1); 0];
        end
      d = 1;
    else
      d = 0;
%       [str2num(File.NT(a).Number) str2num(File.NT(b).Number)]
    end
  else
    if strcmp(Node(n).type,'Initial')
        Node(n).lpar=Node(n).lpar+1;
    else
        Node(n).lpar=Node(n).lpar+[ones(16,1); 0];
    end
    d = 0;
  end

  % ---------------------------------- Show additional interactions a has
 
  if (((d == 2) & (length(Interact{a}.Categ) > 1)) | (d == 1)),
      if BP==1
          n=n+1;
          t=1;
          Node(n).type = 'Motif';
          BP=0;
          Node(n).lpar        = [0.001*ones(16,1); 0];     % left insertion parameters
          Node(n).rpar        = [0.001*ones(16,1); 0];     % right insertion parameters      
          Node(n).IBases(t,:) = [str2num(File.NT(a).Number) str2num(File.NT(b).Number)];
          Node(n).Edge(t)    = c;
      end
    if ~strcmp(Node(n).type,'Initial')
        Node(n).type = 'Motif';
    end
    for k = d:length(Interact{a}.Categ),
      bb = Interact{a}.Index(k);
      cc = Interact{a}.Categ(k);
      t=t+1;
      Node(n).IBases(t,:) = [str2num(File.NT(a).Number) str2num(File.NT(bb).Number)];
      Node(n).Edge(t)    = cc;
    end
  end

  % ---------------------------------- Show additional interactions b has

  if d == 2,                              % if there was a primary interaction
    if length(Interact{b}.Categ) > 1,     % and b makes more than one interact
      if BP==1
          n=n+1;
          t=1;
          Node(n).type = 'Motif';
          BP=0;
          Node(n).lpar        = [0.001*ones(16,1); 0];     % left insertion parameters
          Node(n).rpar        = [0.001*ones(16,1); 0];     % right insertion parameters
          Node(n).IBases(t,:) = [str2num(File.NT(a).Number) str2num(File.NT(b).Number)];
          Node(n).Edge(t)    = c;
      end
      if ~strcmp(Node(n).type,'Initial')
          Node(n).type = 'Motif';
      end
      for k=1:length(Interact{b}.Categ),
        bb = Interact{b}.Index(k);
        cc = Interact{b}.Categ(k);
        if bb ~= a,
          t=t+1;
          Node(n).IBases(t,:) = [str2num(File.NT(b).Number) str2num(File.NT(bb).Number)];
          Node(n).Edge(t)    = cc;           
        end
      end
    end
  end
  a = a + 1;

end


