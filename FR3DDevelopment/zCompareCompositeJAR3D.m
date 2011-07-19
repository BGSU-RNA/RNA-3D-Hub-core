
% File = zAddNTData({'2avy','1j5e'});
if isempty(File(1).Distance),
  for f = 1:length(File),
    c = cat(1,File(f).NT(1:File(f).NumNT).Center); % nucleotide centers
    File(f).Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
  end
end

Verbose = 2;
a = 0;

% -------------------------------------- Load JesseRyan composite alignment

load JesseRyan16SComposite

a = a+1;
Al(a).ModelStructure = JRComposite2AVY;
Al(a).InferStructure = JRComposite1J5E;
Al(a).Name           = 'Jesse-Ryan composite';

% ------------------------------ JAR3D alignment of 1j5e seq to 2avy structure

%>E. coli sequence from 2avy -1454.2420647192616 

ES = 'UGAA---GAGU<UUGAUCAUG>GCUC--AGAUUGAACGCUG--GCGGCAGG-{CCUAACA{CAUGCAA--GUC{GAACG{GUAA{CAGGAAGAAGC----------<UU>----------G-CUUCUUUGCUG}ACG}AGU}GGCGGACG-GGUGAGU-AAU-GUCUGGGAA-ACUGCCUGAUGGA--GGGGGAUAACUACUGG<AA>ACGGUAGCUAAUACCGC--AUAACGUC----------G<CA>A----------GACCAA-A-GAGGGGGACC----------U<UC>G----------GGCCUCUU-GCCAUCGGAUGUGCCCAGAU--GGG{AUUAGCUAGUAGGUGGG<GUAACGG>CUCACCUAGGCGACGAU}CCC-UA-GCUGGUCUG<AG>AGGAUGACCAGC-CACA----CUGGAACU<GAG>ACACGGUCCAG--ACUCC----------U<AC>G----------GGAGG--CAGCAG}UGGGGAAU}--AUUGC{ACAAUGGGC----------G<CA>A----------GCCUGAU}GCAG-CCAUGCCGCG-U-GUAUG{AAGAAGGCC----------U<UC>G----------GGUUGUAA}AGUAC--{UUUCAGCGGGGAGGAAGGGAGUAAAG<UUAAUAC>CUUUGCUC-AUUGACGUUACCCGCAGAAG}-AA---GC{ACCGGCUAACUCCG{UGCCAGC<AGCC>GCGGUAAUA}CGGAG-GGU}GC-AAGCGUUAAUCGG-AAUU--ACUG{GGCGUAAAGC--GCACGCAGGCG-GUUUGUU{AAG{UCAGAUGUGAAAUCCCCGGGC<UCA>ACCUGGGAACUGCAUCUGAUA}CU}GGCAAGC--{UUGAGUCUCGUAGAGGGGGGU-AGAAUUCCAGGUGUAGCG----------G<UGAA>A----------UGCGUAGAGAUCUGGAGGAAUACCG-GUGGC<GAAG>GCG-GCCCCCUGGACGAAGACUGAC}-GCUCAGGUGCG-AAAGCGUGGGGAGCAAACAGG----------<AUUAGAUAC>----------CCUGGUAGUCCACGCCGU--AAACGA{UGUCGACU-UGGAGGUUGUGCC------------<CUUGA>------------GGCGUGGCUUCCGG-AGC<UAAC>GCG-UUAAGUCGAC}-CGCC}U--GGGGA{GUACGGCC----------G<CA>A----------GGUUAA}AACUC-AAAUGAAU----{UG-ACGGGGGCCCGC{ACAA-GCGGUG-GAGCAUGUGGUUUAA-UUCG--------A<UGCA>A--------CGCGAAGAACCUUA-CCUGGUCUU-{GACAUCCACGGAA--GUUUUCA<GAGA>UGAGAAUGU-GCC----------<UUCG>----------GGA--ACCGUGAGAC}--{AGGUGCU{GCAUG{GCUGUCGUCA-GCUCGU-GUU<GUGA>AAUGU-UGGGU----------<UAAG>----------UCCCG-CAACGAGCGCAAC-CCUUA-UCCUUU{GUUGC{CAG{CGG----------<UCCG>----------GCCG}-G}G-AACUC}AAAGGA--{GACUGCCAGUG<AUA>A-ACUGGA-GGAAG}-GUGGGG--AUGACGUCAAGUCAU}CAU}GGC-CCUUAC}-GACCAGG-GCUACACACGUGCUA--C-{AAUGGCGCAUACAAAGAGAA{GCGACCUC<GCGA>GAGCAAGCG}GACCUCAUAAAGUGCGUC-GUAGU}C-CGGAUUGGAGUC----------U<GCA>A----------CUCGACUCCAUGAAG-UCG-GAAUCGCU-A{GUAAUCGUGGAU<CAGA>-AUGCCACGGUGA}A-UACGUUCC}CGGGCCUUGUACACAC}--CGCCCGUCACACCAUGGGAGUGGGUUGCAAAAGAAGUAGG{UAGCUUAA-CC----------U<UC>G----------GGAGGGCG}CUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAG-GU-AACCGUAGG<GGAA>CCUGCGGUU--------GGAUCA';

%>T.th. sequence from 1j5e -1712.0964047104503 
TS = 'UGGA---GAGU<UUGAUCCUG>GCUC--AGGGUGAACGCUG--GCGGCGUG-{CCUAAGA{CAUGCAA--GUC{GUGCG{G---{GCCGCGGGGU-----------<UU>-----------UACUCCGUGGUC}-AG}CGG}CGG--ACG-GGUGAGU-AAC-GCGUGGGUGACCUACCCGGAAGA--GGGGGACAACCCGGGG<AA>ACUCGGGCUAAUCCCCC--AUGU-GGACCCGCCCC---<UU>---GGGGUGUGUCC-AA---AGGG----------------C<UU>U--------------GCCC-GCUUCCGGAUGGGCCCGCGU--CCC{AUCAGCUAGUUGGUGGG<GUAAUGG>CCCACCAAGGCGACGAC}GGG-UA-GCCGGUCUG<AG>AGGAUGGCCGGC-CACA----GGGGCACU<GAG>ACACGGGCCCC--ACUCC----------U<AC>G----------GGAGG--CAGCAG}UUAGGAAU}--CUUCC{GCAAUGGGC----------G<CA>A----------GCCUGAC}GGAG-CGACGCCGCU-U-GGAGG{AAGAAG-CCC---------U<UC>G---------GGG-UGUAA}ACUCC--{UG-AACCCGGGA-C-G-A-AACCC-C<CGACGAG>GGG-ACU-GACG--G--UACCGGGGUAAU}--A---GC{GCCGGCCAACUCCG{UGCCAGC<AGCC>GCGGUAAUA}CGGAG-GGC}GC-GAGCGUUACCCGG-AUUC--ACUG{GGCGUAAAGG--GCGUGUAGGCG-GCCUGGG{GCG{UCCCAUGUGAAAGACCACGGC<UCA>ACCGUGGGGGAGCGUGGGAUA}CG}CUCAGGC--{UAGACGGUGGGAGAGGGUGGU-GGAAUUCCCGGAGUAGCG----------G<UGAA>A----------UGCGCAGAUACCGGGAGGAACGCCG-AUGGC<GAAG>GCA-GCCACCUGGUCCACCCGUGAC}-GCUGAGGCGCG-AAAGCGUGGGGAGCAAACCGG----------<AUUAGAUAC>----------CCGGGUAGUCCACGCCCU--AAACGA{UGCGCGCU-AGGU-CUC-UGG-------------<GUCUC>-------------CUG-GGG-GCCGA-AGC<UAAC>GCG-UUAAGCGCGC}-CGCC}U--GGGGA{GUACGGCC----------G<CA>A----------GGCUGA}AACUC-AAAGGAAU----{UG-ACGGGGGCCCGC{ACAA-GCGGUG-GAGCAUGUGGUUUAA-UUCG--------A<AGCA>A--------CGCGAAGAACCUUA-CCAGGCCUU-{GACAUGCUAGGGAA-CCCGGGU<GAAA>GCCUGGGGU-GCCCC--------<GCGA>--------GGGGA-GCCCUAG-CAC}--{AGGUGCU{GCAUG{GCCGUCGUCA-GCUCGU-GCC<GUGA>GGUGU-UGGGU----------<UAAG>----------UCCCG-CAACGAGCGCAAC-CCCCG-CCGUUA{GUUGC{CAG{CGG----------<UUCG>----------GCCG}-G}GCA-CUC}UAACGG--{GACUGCCCGC-<GAA>-AGCGGGA-GGAAG}--GAGGG-GACGACGUCUGGUCAG}CAU}GGC-CCUUAC}-GGCCUGG-GCGACACACGUGCUA--C-{AAUGCCCACUACAAAGCGAU{GCCACCCG<GCAA>CGGGGAGCU}AAUCGCAAAAAGGUGGGCCCAGUU}--CGGAUUGGGGUC----------U<GCA>A----------CCCGACCCCAUGAAG-CCG-GAAUCGCU-A{GUAAUCGCGGAU<CAGC>CAUGCCGCGGUGA}A-UACGUUCC}CGGGCCUUGUACACAC}--CGCCCGUCACGCCAUGGGAGCGGGCUCUACCCGA----AG{UCGCCGGGAGCC---------U<AC>G---------GGC-AGGCG}CC---GAGGGUAGGGCCCGUGACUGGGGCGAAGUCGUAACAAG-GU-AGCUGUACC<GGAA>GGUGCGGCU--GGAUCACUUUCU';


[i1,i2] = pGetAlignedIndices(ES,TS);

a = a+1;
Al(a).ModelStructure = i1;
Al(a).InferStructure = i2;
Al(a).Name           = 'JAR3D 2002 IM Extensible';

%------------------------------------------------------------------------

% ---------------------------------------- Visually compare to Alignment 1

if Verbose > 1,

for a = 2:length(Al),
  clf
  zCompareAlignment(File,Al(1).ModelStructure,Al(1).InferStructure,Al(a).ModelStructure,Al(a).InferStructure,Al(1).Name,Al(a).Name);
  Titl = ['Agreement between ' Al(a).Name ' and ' Al(1).Name];
%  title(Titl);
  saveas(gcf,[strrep(Titl,' ','_') '.pdf'],'pdf');
end

% --------------------------------------- Calculations for each alignment

for a = 1:length(Al),
  Al(a).Matrix = sparse(Al(a).ModelStructure,Al(a).InferStructure,ones(1,length(Al(a).ModelStructure)));
  Al(a).Matrix(length(File(1).NT),length(File(2).NT)) = 0;

  Comp = File(2);
  Comp.Filename = [Comp.Filename ' inferred from ' NewFile(1).Filename];
  Comp.Edge  = Comp.Edge .* (fix(abs(NewFile(1).Edge)) == fix(abs(NewFile(2).Edge)));
  Comp.BasePhosphate  = Comp.BasePhosphate .* (fix(abs(NewFile(1).BasePhosphate/100)) == fix(abs(NewFile(2).BasePhosphate/100)));

  if Verbose > 2,
    Al(a).Tally = zAlignmentDiagram(File,Al(a).ModelStructure,Al(a).InferStructure,1);
  else
    Al(a).Tally = zAlignmentDiagram(File,Al(a).ModelStructure,Al(a).InferStructure,0);
  end
end

a = 2;
zAlignmentDiagram(File,Al(a).ModelStructure,Al(a).InferStructure,1);


% --------------------------------------- Superimpose local neighborhoods

for a = 1:length(Al),
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

for a = 1:length(Al),
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

for a = 1:length(Al),
  Agree  = sum(sum(Al(a).Matrix .* Al(1).Matrix == 1));
  Missed = sum(sum(Al(1).Matrix > Al(a).Matrix));
  Extra  = sum(sum(Al(a).Matrix > Al(1).Matrix));
  Identical = length(find(cat(1,File(1).NT(Al(a).ModelStructure).Code)==cat(1,File(2).NT(Al(a).InferStructure).Code)));

  T{ 1,a+1} = Al(a).Name;
  T{ 2,a+1} = length(Al(a).ModelStructure);
  for v = 1:6,
    T{v+2,a+1} = Al(a).Tally(v);
  end
  T{ 9,a+1} = Agree;
  T{10,a+1} = Missed;
  T{11,a+1} = Extra;
  T{12,a+1} = mean(Al(a).Discrep);
  T{13,a+1} = median(Al(a).Discrep);
  T{14,a+1} = mean(Al(a).IDI);
  T{15,a+1} = median(Al(a).IDI);
  T{16,a+1} = Identical;
end

T{1,1} = [File(1).Filename ' as the model, ' File(2).Filename ' as the unknown structure'];
T{2,1} = 'Number aligned';
T{3,1} = ['Nested cWW correctly inferred ' File(2).Filename];
T{4,1} = ['Nested non-cWW correctly inferred in ' File(2).Filename];
T{5,1} = ['Non-nested cWW correctly inferred in ' File(2).Filename];
T{6,1} = ['Non-nested non-cWW correctly inferred in ' File(2).Filename];
T{7,1} = ['Stacking correctly inferred in ' File(2).Filename];
T{8,1} = ['Base-phosphate correctly inferred in ' File(2).Filename];
T{9,1} = ['Number of correspondences agreeing with ' Al(1).Name];
T{10,1} = ['Number of correspondences missing, compared to ' Al(1).Name];
T{11,1} = ['Number of correspondences extra, compared to ' Al(1).Name];
T{12,1} = 'Mean geometric discrepancy at 8 Angstroms';
T{13,1} = 'Median geometric discrepancy at 8 Angstroms';
T{14,1} = 'Mean IDI between aligned base combinations and real pairs';
T{15,1} = 'Median IDI between aligned base combinations and real pairs';
T{16,1} = 'Number of exact base matches in alignment';
xlswrite('16S_alignment_comparison.xls',T);

T


end
