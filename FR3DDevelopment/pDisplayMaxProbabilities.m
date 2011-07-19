% pDisplayMaxProbabilities(Node,Sequence) gives a graphical display of 
% maximum probabilities for each node

function [void] = pDisplayMaxProbabilities(Node,Sequence)

K = length(Sequence);
N = length(Node);

for k=1:K,                              % loop through sequences
 for n=1:N,                               % loop through nodes
  P = Sequence(k).P;

  if ~isempty(P),                        % if sequence has been parsed
    clf
    L  = length(Sequence(k).X);
    MP = -Inf*ones(L,L);                 % matrix of -infinity
 
    for s=1:Node(n).numstates,           % find state w/ highest probability
      MP = max(MP,P(n,s).mp);
    end

    MP = MP + tril(-Inf*ones(L,L));      % make lower triangular part -Inf    
           
    subplot(2,1,1);

%    pcolor(MP);
    mesh(1:L,1:L,MP');      
    view(2);

    hold on

    iMin = Sequence(k).Constraint(n).iMin;
    iMax = Sequence(k).Constraint(n).iMax;
    jMin = Sequence(k).Constraint(n).jMin;
    jMax = Sequence(k).Constraint(n).jMax;

    plot([iMin iMax iMax iMin iMin],[jMin jMin jMax jMax jMin]);

    i = Sequence(k).TraceInfo(n).i;
    j = Sequence(k).TraceInfo(n).j;

    plot(i,j,'k.','markersize',18);

    fprintf('Node %d generated subsequence %d to %d\n', n, i, j);

    title([strrep(Sequence(k).Organism,'_','\_') ' Node ' num2str(n)]);
    xlabel('Starting index of subsequence');
    ylabel('Ending index');

    for d = 1:(L-1),
      MPM(d) = max(diag(MP,d));
    end

    minval = min(MPM(find(MPM > -Inf)));

    subplot(2,1,2)
    plot(MPM-max(MPM))
    grid on
    axis([0 L minval 0]);

    title([strrep(Sequence(k).Organism,'_','\_') ' Node ' num2str(n)]);
    xlabel('Max over subsequences with specified length, versus length');

    fprintf('Node %d type %s\n', n, Node(n).type);

    pause
  end
 end
end

