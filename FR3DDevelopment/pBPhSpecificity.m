% pBPhSpecificity(bph,bp) returns a 4x1 probability distribution over ACGU telling how likely each of these is, given that the nucleotide makes BPh interaction with FR3D code bph.  If the nucleotide is part of a basepair, set bp = 1; some day we can calculate Q differently in this case

function [Q] = pBPhSpecificity(bph,bp)

if nargin < 2,
  bp = 1;
end

Code = [1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 2 3]; % which base makes which BPh

Q = 0.10*ones(4,1);
Q(Code(bph),1) = 0.7;                   % 70% conservation of the base
                                        % is about the lowest observed in any
                                        % combination of circumstances

if bp == 1,                             % the base is part of a pair

  % more specific calculations can go here

else


end

