
% ===========================================================================
% ===========================================================================
% ===========================================================================
% ===========================================================================
% ===========================================================================
% ===========================================================================
% ===========================================================================

% old code below

return




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
                Node(n).type = 'Cluster';
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
  

% probe for insertions

    if strcmp(Node(n).type,'Initial')
      Node(n).lpar=Node(n).lpar+1;
    else
      Node(n).lpar=Node(n).lpar+[ones(16,1); 0];
    end
    d = 0;

  % ---------------------------------- Show additional interactions a has
 
  if (((d == 2) & (length(Interact{a}.Categ) > 1)) | (d == 1)),
      if BP==1
          n=n+1;
          t=1;
          Node(n).type = 'Cluster';
          BP=0;
          Node(n).lpar        = [0.001*ones(16,1); 0];     % left insertion parameters
          Node(n).rpar        = [0.001*ones(16,1); 0];     % right insertion parameters      
          Node(n).IBases(t,:) = [str2num(File.NT(a).Number) str2num(File.NT(b).Number)];
          Node(n).Edge(t)    = c;
      end
    if ~strcmp(Node(n).type,'Initial')
        Node(n).type = 'Cluster';
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
          Node(n).type = 'Cluster';
          BP=0;
          Node(n).lpar        = [0.001*ones(16,1); 0];     % left insertion parameters
          Node(n).rpar        = [0.001*ones(16,1); 0];     % right insertion parameters
          Node(n).IBases(t,:) = [str2num(File.NT(a).Number) str2num(File.NT(b).Number)];
          Node(n).Edge(t)    = c;
      end
      if ~strcmp(Node(n).type,'Initial')
          Node(n).type = 'Cluster';
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
