
AnalyzeList = [1:8];
AnalyzeList = [1:11];

File = zAddNTData({'2avy','1j5e'});
if isempty(File(1).Distance),
  for f = 1:length(File),
    c = cat(1,File(f).NT(1:File(f).NumNT).Center); % nucleotide centers
    File(f).Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
  end
end

Verbose = 1;
a = 0;

% -------------------------------------- Load JesseRyan composite alignment

load JesseRyan16SComposite

a = a+1;
Al(a).ModelStructure = JRComposite2AVY;
Al(a).InferStructure = JRComposite1J5E;
Al(a).Name           = 'Jesse-Ryan composite';

% -------------------------------------- Load Ryan's 3D to 3D alignment

load 16SalignmentInd

a = a+1;
Al(a).ModelStructure = ind2';
Al(a).InferStructure = ind1';
Al(a).Name           = 'Ryan 3D to 3D';

% -------------------------------------- Load Jesse's 3D to 3D alignment

[File1,i1,File2,i2] = zReadJesseAlignments('16S');

a = a+1;
Al(a).ModelStructure = i1';
Al(a).InferStructure = i2';
Al(a).Name           = 'Jesse 3D to 3D';

% ------------------------------ JAR3D1 alignment of 1j5e seq to 2avy structure

load JAR3D1_Alignment

a = a+1;
Al(a).ModelStructure = i1;
Al(a).InferStructure = i2;
Al(a).Name           = 'JAR3D1';

% ------------------------------ JAR3D2 alignment of 1j5e seq to 2avy structure

load JAR3D2_Alignment

a = a+1;
Al(a).ModelStructure = i1;
Al(a).InferStructure = i2;
Al(a).Name           = 'JAR3D2';

% ------------------------------ JAR3D3 alignment of 1j5e seq to 2avy structure

load JAR3D3_Alignment

a = a+1;
Al(a).ModelStructure = i1;
Al(a).InferStructure = i2;
Al(a).Name           = 'JAR3D3';

% ------------------------------ JAR3D4 alignment of 1j5e seq to 2avy structure

load JAR3D4_Alignment

a = a+1;
Al(a).ModelStructure = i1;
Al(a).InferStructure = i2;
Al(a).Name           = 'JAR3D4';

% ------------------------------ JAR3D5 alignment of 1j5e seq to 2avy structure

load JAR3D5_Alignment

a = a+1;
Al(a).ModelStructure = i1;
Al(a).InferStructure = i2;
Al(a).Name           = 'JAR3D5';

% ------------------------------ JAR3D6 alignment of 1j5e seq to 2avy structure

load JAR3D6_Alignment

a = a+1;
Al(a).ModelStructure = i1;
Al(a).InferStructure = i2;
Al(a).Name           = 'JAR3D6';

% --------------------------------------- Needleman-Wunsch alignment

for w = 1:2,
  if any(a+1 == AnalyzeList),

%    [matches,align1,align2,s1,s2] = zNeedlemanWunsch(cat(2,File(1).NT.Base),cat(2,File(2).NT.Base),2*w);
    g = -6.56 + 0.02*w;
    [matches,align1,align2,s1,s2] = dNeedlemanWunsch(cat(2,File(1).NT.Base),cat(2,File(2).NT.Base),0.8,g);      % Ryan's alignment has 80% identical

    a = a+1;
    Al(a).ModelStructure = align1;
    Al(a).InferStructure = align2;
    Al(a).Name           = ['NW del ' num2str(g)];
  end
end

%------------------------------------------------------------------------
% --------------------------------------- Calculations for each alignment

for a = AnalyzeList,
  Al(a).Matrix = sparse(Al(a).ModelStructure,Al(a).InferStructure,ones(1,length(Al(a).ModelStructure)));

  Al(a).Matrix(length(File(1).NT),1) = 0;    % make the matrix bigger
  Al(a).Matrix(1,length(File(2).NT)) = 0;

  if Verbose > 1,
    Al(a).Tally = zAlignmentDiagram(File,Al(a).ModelStructure,Al(a).InferStructure,1);
  else
    Al(a).Tally = zAlignmentDiagram(File,Al(a).ModelStructure,Al(a).InferStructure,0);
  end
end

% --------------------------------------- Superimpose local neighborhoods

for a = AnalyzeList,
  Al(a).Discrep = zHistogramDiscrepanciesInAlignment(File,Al(a).ModelStructure,Al(a).InferStructure);
end

figure(10)
clf

m = 3;
s = ceil(sqrt(length(Al)));
for a = 1:length(Al),
  subplot(s,s,a)
  if length(Al(a).Discrep) > 0,
    n = hist(min(Al(a).Discrep,m),-0.025+(0:0.05:m));
    hist(min(Al(a).Discrep,m),-0.025+(0:0.05:m))
    axis([0 m 0 max(n)*1.1])
  end
  title(Al(a).Name);
end

saveas(gcf,'Histogram Discrepancies in 16S Alignments.pdf','pdf')

% -------------------------------- Calculate IDI of aligned base combinations 
%                                  with basepairs in 3D structure

for a = AnalyzeList,
  Al(a).IDI     = zAlignmentIsostericity(File,Al(a).ModelStructure,Al(a).InferStructure);
end

figure(11)
clf

m = 10;
s = ceil(sqrt(length(Al)));
for a = 1:length(Al),
  subplot(s,s,a)
  if length(Al(a).IDI) > 0,
    n = hist(min(Al(a).IDI,m),-0.025+(0:0.05:m));
    hist(min(Al(a).IDI,m),-0.025+(0:0.05:m))
    axis([-1 m 0 max(n)*1.1])
  end
  title(Al(a).Name);
end

saveas(gcf,'Histogram of IDI values in 16S Alignments.pdf','pdf')




% --------------------------------------- Compare to Alignment 1

for a = AnalyzeList,
%for a = 1:1,
  Agree  = sum(sum(Al(a).Matrix .* Al(1).Matrix == 1));  % where they agree
  Missed = sum(sum(Al(1).Matrix > Al(a).Matrix));
  Extra  = sum(sum(Al(a).Matrix > Al(1).Matrix));
  Identical = length(find(cat(1,File(1).NT(Al(a).ModelStructure).Code)==cat(1,File(2).NT(Al(a).InferStructure).Code)));


  % ----------------- Find identical pairs that also correspond to alignment 1

  [i,j,k] = find(Al(a).Matrix .* Al(1).Matrix); % part of alignment that agrees

  Al1 = i';
  Al2 = j';

  E1 = fix(File(1).Edge);                         % interactions in File 1
  E1 = E1 - 2*E1.*((E1 == -1) + (E1 == -2) + (E1 == -7) + (E1 == -8) + (E1 == -11) + (E1 == -12));                            % use abs values for these pairs

  E2 = fix(File(2).Edge);
  E2 = E2 - 2*E2.*((E2 == -1) + (E2 == -2) + (E2 == -7) + (E2 == -8) + (E2 == -11) + (E2 == -12));

  FF = File(2);
  FF.Edge = FF.Edge(Al2,Al2) .* (E1(Al1,Al1) == E2(Al2,Al2));

%  Tally = zSummarizeInteractions(FF,0);

  Tally = zAlignmentDiagram(File,i',j',Verbose);

  T{ 1,a+1} = Al(a).Name;
  T{ 2,a+1} = length(Al(a).ModelStructure);
  for v = 1:6,
    T{v+2,a+1} = Tally(v);
  end
  T{ 9,a+1} = Agree;
  T{10,a+1} = Missed;
  T{11,a+1} = Extra;
  T{12,a+1} = mean(Al(a).Discrep);
  T{13,a+1} = median(Al(a).Discrep);
  T{14,a+1} = mean(Al(a).IDI);
  T{15,a+1} = median(Al(a).IDI);
  T{16,a+1} = Identical;
  for v = 1:6,
    T{v+16,a+1} = Al(a).Tally(v);
  end
end

T{ 1,1} = [File(1).Filename ' as the model, ' File(2).Filename ' as the unknown structure'];
T{ 2,1} = 'Number aligned';
T{ 3,1} = ['Nested cWW correctly inferred and placed in' File(2).Filename];
T{ 4,1} = ['Nested non-cWW correctly inferred and placed in ' File(2).Filename];
T{ 5,1} = ['Non-nested cWW correctly inferred and placed in ' File(2).Filename];
T{ 6,1} = ['Non-nested non-cWW correctly inferred and placed in ' File(2).Filename];
T{ 7,1} = ['Stacking correctly inferred and placed in ' File(2).Filename];
T{ 8,1} = ['Base-phosphate correctly inferred and placed in ' File(2).Filename];
T{ 9,1} = ['Number of correspondences agreeing with ' Al(1).Name];
T{10,1} = ['Number of correspondences missing, compared to ' Al(1).Name];
T{11,1} = ['Number of correspondences extra, compared to ' Al(1).Name];
T{12,1} = 'Mean geometric discrepancy at 8 Angstroms';
T{13,1} = 'Median geometric discrepancy at 8 Angstroms';
T{14,1} = 'Mean IDI between aligned base combinations and real pairs';
T{15,1} = 'Median IDI between aligned base combinations and real pairs';
T{16,1} = 'Number of exact base matches in alignment';
T{17,1} = 'Nested cWW correctly inferred, maybe not placed right ';
T{18,1} = 'Nested non-cWW correctly inferred, maybe not placed right ';
T{19,1} = 'Non-nested cWW correctly inferred, maybe not placed right ';
T{20,1} = 'Non-nested non-cWW correctly inferred, maybe not placed right ';
T{21,1} = 'Stacking correctly inferred, maybe not placed right ';
T{22,1} = 'Base-phosphate correctly inferred, maybe not placed right ';

xlswrite('16S_alignment_comparison.xls',T);

T

% ---------------------------------------- Visually compare to Alignment 1

if Verbose > 0,

  for a = AnalyzeList,
    if a ~= 1,
      clf
      zCompareAlignment(File,Al(1).ModelStructure,Al(1).InferStructure,Al(a).ModelStructure,Al(a).InferStructure,Al(1).Name,Al(a).Name);
      Titl = ['Agreement between ' Al(a).Name ' and ' Al(1).Name];
%  title(Titl);
      saveas(gcf,[strrep(Titl,' ','_') '.pdf'],'pdf');
    end
end

end
