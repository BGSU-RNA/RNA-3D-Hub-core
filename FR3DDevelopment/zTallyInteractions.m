% zTallyInteractions(File) counts the number of interactions in the files in the variable file, adding each of them as a row of the variable Tally.  zTallyInteractions(File,A1,A2) calculates the number of conserved interactions using aligned nucleotides A1 from File(1) and A2 from File(2)

function [Tally] = zTallyInteractions(File,A1,A2,Verbose)

if nargin == 1,                         % single file

  for f = 1:length(File),

    E  = fix(abs(File(f).Edge));
    B  = E .* (E > 0) .* (E < 24);        % pairs and stacks
    C  = File(f).Crossing;

    Color = (B==1).*(C==0) + 2*(B>1).*(B<13).*(C==0) + 3*(B==1).*(C>0) + 4*(B > 1).*(B < 13) .*(C>0) + 5*(B > 20) .* (B < 25);
                                      % disjoint possibilities

    BP = abs(File(f).BasePhosphate);         % BPh interactions
    BP = (BP > 0) .* (BP < 100);          % exclude near BPh and self interactions

    [i,j,c] = find(triu(Color));          % pull out indices and interactions
 
    [ii,jj,cc] = find(6*BP);              % handle BPh separately, don't add

    k = find(ii ~= jj);                   % eliminate self interactions

    i = [i; ii(k)];                       % append BPh to other interactions
    j = [j; jj(k)];
    c = [c; cc(k)];

    cww = length(find(c == 1));
    Tally(f,1) = cww;
    noncww = length(find(c == 2));
    Tally(f,2) = noncww;
    nonnestcww = length(find(c == 3));
    Tally(f,3) = nonnestcww;
    nonnestnoncww = length(find(c == 4));
    Tally(f,4) = nonnestnoncww;
    stack = length(find(c == 5));
    Tally(f,5) = stack;
    bph = length(find(c == 6));
    Tally(f,6) = bph;
  end
else                                    % tally conserved interactions

  E1 = fix(File(1).Edge(A1,A1));         % pairwise interactions among aligned

  E1 = E1 - 2*E1.*((E1 == -1) + (E1 == -2) + (E1 == -7) + (E1 == -8) + (E1 == -11) + (E1 == -12));                            % use abs values for these pairs

  E2 = fix(File(2).Edge(A2,A2));
  E2 = E2 - 2*E2.*((E2 == -1) + (E2 == -2) + (E2 == -7) + (E2 == -8) + (E2 == -11) + (E2 == -12));

  BP1 = fix(File(1).BasePhosphate(A1,A1));
  BP2 = fix(File(2).BasePhosphate(A2,A2));

  BP1 = BP1 .* (abs(BP1) > 0) .* (abs(BP1) < 100);  % only use true BPh
  BP2 = BP2 .* (abs(BP2) > 0) .* (abs(BP2) < 100);  % only use true BPh

  NewFile.Edge = E1 .* (E1 == E2);          % conserved interactions
  NewFile.BasePhosphate = BP1 .* (BP1 == BP2);  % conserved BPh


  NewFile.Crossing = File(1).Crossing(A1,A1);
  NewFile.Crossing = File(2).Crossing(A2,A2);

  Tally = zTallyInteractions(NewFile);  
end
