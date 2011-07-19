% pMakeNodesCheckForExtensibility checks the rest of the stem to see if there are long-range basepairs, BPh, or GU packing interactions

Outside = [1:(a-1) (B+1):length(File.NT)];

if sum(sum(J(a:B,Outside))) == 0 && ... % no tertiary inter outside this stem
   sum(sum(File.BasePhosphate(a:B,Outside)>0)) == 0 && ... % no BPh either
   sum(HasGUPacking(a:B)) == 0 && ...
   (~strcmp(File.Filename,'2J01') || (a > 2842)) && ...
   (~strcmp(File.Filename,'2AW4') || (a > 69)),

  if TertiaryFreeNode == 0,              % only mark the first one
    TertiaryFreeNode = n;                
    if Verbose > 0,
     fprintf('After node %3d, this stem is free of tertiary interactions\n',n);
    end

    Node(1).Extensibility(a:B) = ones(1,length(a:B));

    if Extension > 1,                    % replace rest of stem with cWW BPs
      G(a:B,a:B) = G(a:B,a:B) .* (abs(fix(G(a:B,a:B))) == 1);
      H(a:B,a:B) = H(a:B,a:B) .* (abs(fix(G(a:B,a:B))) == 1);
      HasMotif(a:B) = zeros(1,length(a:B));
    end
  end
end

