% pLCS.m finds longest common substrings of x in the strings in array Y,
% Only substrings of length >= minlength are returned.  
% Substrings that are nested within larger substrings are not returned.

function [FinalMatch] = zLCS1(x,Y,minlength)

Lx = length(x);

n = 1;

for s = 1:length(Y),
  y  = Y{s};
  Ly = length(y);
  
  a = [Lx:-1:1 ones(1,Ly-1)];
  b = [Lx*ones(1,Ly) (Lx-1):-1:1];
  c = [ones(1,Lx) 2:Ly];
  d = [1:Ly Ly*ones(1,Lx-1)];
  
  h = abs(a-c) + abs(b-d);            
  [m, i] = sort(h);                  % start with most overlapping

  for p = 1:5,
    S = x(a(i(p)):b(i(p)));
    T = y(c(i(p)):d(i(p)));
    g = [0 (S==T) 0];               % locations of matches
    h = diff(g);                    % starts and stops of matches
    j = find(h == 1);               % locations of starts
    k = find(h == -1)-1;            % locations of stops
  
    for r = 1:length(j),
      if k(r)-j(r)+1 >= minlength,    
        Match(n).Length = k(r)-j(r)+1;
        Match(n).XStart  = j(r)+a(i(p))-1;
        Match(n).XStop   = k(r)+a(i(p))-1;
        Match(n).YStart  = j(r)+c(i(p))-1;
        Match(n).YStop   = k(r)+c(i(p))-1;
        Match(n).Seqnum = s;
        n = n + 1;
      end
    end
  end
end

if exist('Match') > 0,
  lengths = cat(1,Match(:).Length);
  [m, i] = sort(-lengths'+0.0001*(1:length(lengths)));
  Match = Match(i);

  FinalMatch(1) = Match(1);
  f = 1;
  
  for r = 2:length(Match),
    M = Match(r);
%    fprintf('Seq %2d Length %2d %2d to %2d %s\n', M.Seqnum, M.Length, M.XStart, M.XStop, x(M.XStart:M.XStop));
  
    Omit = 0;
    for g = 1:f,
     FM = FinalMatch(g);
      if ((Match(r).XStart >= FM.XStart) & (Match(r).XStop <= FM.XStop)),
        Omit = 1;
      end
    end

    if Omit == 0, 
      f = f + 1;
      FinalMatch(f) = Match(r);
    end
    
  end
else
  FinalMatch = [];
end

%fprintf('\n');

for r = 1:length(FinalMatch),
  M = FinalMatch(r);
%  fprintf('Seq %2d Length %2d %2d to %2d %s\n', M.Seqnum, M.Length, M.XStart, M.XStop, x(M.XStart:M.XStop));
end
