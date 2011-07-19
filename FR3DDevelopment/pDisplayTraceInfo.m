% pDisplayTraceInfo(Node,Sequence) gives a graphical display of nodes
% and the subsequences they explain

function [void] = pDisplayTraceInfo(Node,Sequence)

K = length(Sequence);
N = length(Node);

for k=1:K,
  t = Sequence(k).TraceInfo;

  if ~isempty(t),
    clf
    plot(t(1).i,t(1).j,'*');
    hold on
    plot([t(1).i t(1).j], [t(1).i t(1).j]);

    for n=2:N,
      plot(t(n).i,t(n).j,'*');
      switch Node(n).type
        case {'Initial','Generic','Junction'}, 
          for a = 1:length(Node(n).nextnode),
            nn = Node(n).nextnode(a);
            plot([t(n).i t(nn).i], [t(n).j t(nn).j]);
          end
      end
    end

    axis([0 t(1).j 0 t(1).j]);

    title(strrep(Sequence(k).Organism,'_','\_'));

    pause
  end
end

