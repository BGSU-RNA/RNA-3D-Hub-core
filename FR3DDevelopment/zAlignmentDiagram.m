% zAlignmentDiagram takes an alignment of the NTs in File(1) and File(2) and produces three circular diagrams.  The first two show the interactions in File(1) and File(2), with nucleotides spaced around the circle so that nucleotides in the same location are aligned.  The third circular diagram shows the 

function [Tally] = zAlignmentDiagram(File,Aligned1,Aligned2,Verbose)

if nargin < 4,
  Verbose = 1;
end

% ----------- Make functions from File indices to rows of the alignment
% R1 is a function from the indices of File(1) to the row of the alignment
% Needed to represent the circular diagrams for the two structures in
% the same way, putting aligned nucleotides in the same places around the
% circle

clear R1 R2

m = max(length(File(1).NT),length(File(2).NT));

A1 = [0 Aligned1 m+1];               % add some alignments to simplify
A2 = [0 Aligned2 m+1];

R1(1) = 0;                           % index 1 maps to row 0; fixed below
R2(1) = 0;

for a = 1:(length(A1)-1),                  % run through all correpondences
  r  = max(R1(end),R2(end))+1;             % current row of the alignment
  R1 = [R1 r];                             % these are matched; append to each
  R2 = [R2 r];

  a1 = r + (1:(A1(a+1)-1-A1(a)));          % new rows for 1
  a2 = R1(end)+(1:(A2(a+1)-1-A2(a)));      % new rows for 2

  R1 = [R1 a1];                            % add rows for non-aligned from 1
  R2 = [R2 a2];                            % add rows for non-aligned from 2
end

R1 = R1(3:end)-1;                   % omit first two assignments
R2 = R2(3:end)-1;

R1 = R1(1:length(File(1).NT));      % only keep enough for the actual indices
R2 = R2(1:length(File(2).NT));

% ---------- Make new File variables containing only aligned nucleotides

clear NewFile

m = max(max(R1),max(R2));

N1 = File(1).NT(1);
N1.Number = '';
N1.Base = '';

for f = 1:2,
  NewFile(f).Filename      = File(f).Filename;
  NewFile(f).Edge          = sparse(zeros(m,m));
  NewFile(f).Crossing      = sparse(zeros(m,m));
  NewFile(f).BasePhosphate = sparse(zeros(m,m));
  NewFile(f).Covalent      = sparse(zeros(m,m));
  for n = 1:m,
    NewFile(f).NT(n) = N1;
  end
end

for r = 1:length(R1),
  NewFile(1).NT(R1(r)) = File(1).NT(r);
end

NewFile(1).Edge(R1,R1)          = File(1).Edge;
NewFile(1).Crossing(R1,R1)      = File(1).Crossing;
NewFile(1).BasePhosphate(R1,R1) = File(1).BasePhosphate;
NewFile(1).Covalent(R1,R1)      = File(1).Covalent;

for r = 1:length(R2),
  NewFile(2).NT(R2(r)) = File(2).NT(r);
end

NewFile(2).Edge(R2,R2)          = File(2).Edge;
NewFile(2).Crossing(R2,R2)      = File(2).Crossing;
NewFile(2).BasePhosphate(R2,R2) = File(2).BasePhosphate;
NewFile(2).Covalent(R2,R2)      = File(2).Covalent;

% ------------------------------------ Decide what to display

thickness = 0.1;                     % thickness of arcs in basepair diagram

if Verbose > 0,
  ViewList = [1 1 1 1 1 1 0 0 1 0 0];

  % ------------------------------------ Show interactions in File(1)

  figure(1)
  clf
  zCircularDiagram(NewFile(1),thickness,ViewList);
  saveas(gcf,[NewFile(1).Filename '_diagram.pdf'],'pdf');

  % ------------------------------------ Show interactions in File(1)

  figure(2)
  clf
  zCircularDiagram(NewFile(2),thickness,ViewList);
  saveas(gcf,[NewFile(2).Filename '_diagram.pdf'],'pdf');

  % ------------------------------------ Look for pairs that agree

  E1 = fix(NewFile(1).Edge);                         % interactions in File 1
  E1 = E1 - 2*E1.*((E1 == -1) + (E1 == -2) + (E1 == -7) + (E1 == -8) + (E1 == -11) + (E1 == -12));                            % use abs values for these pairs

  E2 = fix(NewFile(2).Edge);
  E2 = E2 - 2*E2.*((E2 == -1) + (E2 == -2) + (E2 == -7) + (E2 == -8) + (E2 == -11) + (E2 == -12));

  BP1 = fix(NewFile(1).BasePhosphate);
  BP2 = fix(NewFile(2).BasePhosphate);

  % ------------------------------------ Show interactions in File(1) also in 2

  figure(3)
  clf
  Comp = NewFile(1);
%  Comp.Filename = [Comp.Filename ' inferred from ' NewFile(2).Filename];
  Comp.Filename = [Comp.Filename ' conserved'];

  Comp.Edge  = Comp.Edge .* (E1 == E2);

  Comp.BasePhosphate  = Comp.BasePhosphate .* (BP1 == BP2);

  zCircularDiagram(Comp,thickness,ViewList);
  saveas(gcf,[NewFile(1).Filename '-' NewFile(2).Filename '_alignment_conserved.pdf'],'pdf');

  % ------------------------------------ Show interactions in File(2) also in 1

  figure(4)
  clf
  Comp = NewFile(2);
%  Comp.Filename = [Comp.Filename ' inferred from ' NewFile(1).Filename];
  Comp.Filename = [Comp.Filename ' conserved'];

  Comp.Edge  = Comp.Edge .* (E1 == E2);

  Comp.BasePhosphate  = Comp.BasePhosphate .* (BP1 == BP2);

  Tally = zTallyInteractions(Comp);

  zCircularDiagram(Comp,thickness,ViewList);
  saveas(gcf,[NewFile(2).Filename '-' NewFile(1).Filename '_alignment_conserved.pdf'],'pdf');

  % ------------------------------------ Show interactions in File(1) also in 2

  figure(5)
  clf
  Comp = NewFile(1);
  Comp.Filename = [Comp.Filename ' not conserved'];
  Comp.Edge  = Comp.Edge .* (fix(abs(NewFile(1).Edge)) ~= fix(abs(NewFile(2).Edge)));
  Comp.BasePhosphate  = Comp.BasePhosphate .* (fix(abs(NewFile(1).BasePhosphate/100)) ~= fix(abs(NewFile(2).BasePhosphate/100)));
  zCircularDiagram(Comp,thickness,ViewList);
  saveas(gcf,[NewFile(1).Filename '-' NewFile(2).Filename '_alignment_not_conserved.pdf'],'pdf');

  % ------------------------------------ Show interactions in File(2) also in 1

  figure(6)
  clf
  Comp = NewFile(2);
  Comp.Filename = [Comp.Filename ' not conserved'];
  Comp.Edge  = Comp.Edge .* (fix(abs(NewFile(1).Edge)) ~= fix(abs(NewFile(2).Edge)));
  Comp.BasePhosphate  = Comp.BasePhosphate .* (fix(abs(NewFile(1).BasePhosphate/100)) ~= fix(abs(NewFile(2).BasePhosphate/100)));
  zCircularDiagram(Comp,thickness,ViewList);
  saveas(gcf,[NewFile(2).Filename '-' NewFile(1).Filename '_alignment_not_conserved.pdf'],'pdf');


else

%  Tally = zTallyInteractions(File,Aligned1,Aligned2,Verbose);

  Comp = NewFile(2);
  Comp.Filename = [Comp.Filename ' inferred from ' NewFile(1).Filename];
  Comp.Edge  = Comp.Edge .* (fix(abs(NewFile(1).Edge)) == fix(abs(NewFile(2).Edge)));
  Comp.BasePhosphate  = Comp.BasePhosphate .* (fix(abs(NewFile(1).BasePhosphate/100)) == fix(abs(NewFile(2).BasePhosphate/100)));

  Tally = zTallyInteractions(Comp);

end

