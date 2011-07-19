% pSequenceVsStructureAlignment compares Jesse's 3D to 3D alignment to the
% output of the Java parser's alignment

% 1j5e is first, 2avy is second

%             [Alignment,File] = zReadJesse16SAlignment(1,2,'3DAlignment\16S_Ec_Tt_Struct_alignment_8_14_07.xls');

AMatrix1 = sparse(zeros(File(1).NumNT,File(2).NumNT));

for i = 1:length(Alignment.Correspondence),
  for j = 1:length(Alignment.Correspondence(i).File1),
    AMatrix1(Alignment.Correspondence(i).File1(j),Alignment.Correspondence(i).File2(j)) = 1;
  end
end

% The next few lines are specific to Jesse's 1j5e to 2avy alignment of the 16S

[i,j] = find(AMatrix1);

a = find(diff(i) < 0);          % places where something is out of order

for k = 1:length(a),
  ii = i(a(k)+1);
  jj = find(AMatrix1(ii,:));
  for m = 1:length(jj),
    jjj = jj(m);
    AMatrix1(ii,jjj)   = 0;            % wacky alignment pair
    fprintf('Because of an inversion in nucleotide order in %s, nucleotides %s%s in %s and %s%s in %s are unlikely to align\n', File(1).Filename, File(1).NT(ii).Base, File(1).NT(ii).Number, File(1).Filename, File(2).NT(jjj).Base, File(2).NT(jjj).Number, File(2).Filename);
  end
end

b = find(diff(j) < 0);          % places where something is out of order

for k = 1:length(b),
  jj = j(b(k)+1);
  ii = find(AMatrix1(:,jj));
  for m = 1:length(ii),
    iii = ii(m);
    AMatrix1(iii,jj)   = 0;            % wacky alignment pair
    fprintf('Because of an inversion in nucleotide order in %s, nucleotides %s%s in %s and %s%s in %s are unlikely to align\n', File(2).Filename, File(1).NT(iii).Base, File(1).NT(iii).Number, File(1).Filename, File(2).NT(jj).Base, File(2).NT(jj).Number, File(2).Filename);
  end
end

d = find(sum(AMatrix1)>1);          % multiple alignments

for i = 1:length(d),
  e = find(AMatrix1(:,d(i)));
  for j = 2:length(e),
    AMatrix1(e(j),d(i)) = 0;
    fprintf('Nucleotide %s%s in %s has multiple alignments\n', File(2).NT(d(i)).Base, File(2).NT(d(i)).Number, File(2).Filename);
  end
end

d = find(sum(AMatrix1')>1);          % multiple alignments

for i = 1:length(d),
  e = find(AMatrix1(d(i),:));
  for j = 2:length(e),
    AMatrix1(d(i),e(j)) = 0;
    fprintf('Nucleotide %s%s in %s has multiple alignments\n', File(1).NT(d(i)).Base, File(1).NT(d(i)).Number, File(1).Filename);
  end
end

% Display the 3D to 3D alignment horizontally ----------------------------

G1 = '';
G2 = '';
G3 = '';
G4 = '';

[s,t] = size(AMatrix1);

a = 1;
b = 1;

while (a < s) && (b < t),
  if any(abs(File(1).Edge(a,a:min(s,(a+8)))) == 1),
    G3 = [G3 '<'];
  elseif any(abs(File(1).Edge(a,max(1,(a-8)):a)) == 1),
    G3 = [G3 '>'];
  else
    G3 = [G3 ' '];
  end

  if any(abs(File(2).Edge(b,b:min(t,(b+8)))) == 1),
    G4 = [G4 '<'];
  elseif any(abs(File(2).Edge(b,max(1,(b-8)):b)) == 1),
    G4 = [G4 '>'];
  else
    G4 = [G4 ' '];
  end

  if AMatrix1(a,b) == 1,
    G1 = [G1 File(1).NT(a).Base];
    G2 = [G2 File(2).NT(b).Base];
    a  = a + 1;
    b  = b + 1;
  elseif sum(AMatrix1(a,b:end)) > 0,
    G1 = [G1 '-'];
    G2 = [G2 File(2).NT(b).Base];
    b  = b + 1;
  elseif sum(AMatrix1(a:end,b)) > 0,
    G1 = [G1 File(1).NT(a).Base];
    G2 = [G2 '-'];
    a  = a + 1;
  else
    G1 = [G1 File(1).NT(a).Base];
    G2 = [G2 '='];
    G1 = [G1 '='];
    G2 = [G2 File(2).NT(b).Base];
    G3 = [G3 ' '];
    G4 = [G4 ' '];
    a  = a + 1;
    b  = b + 1;
  end

end



fprintf('%s %s.\n', File(1).Filename, G3);
fprintf('%s %s.\n', File(1).Filename, G1);
fprintf('%s %s.\n', File(2).Filename, G2);
fprintf('%s %s.\n', File(2).Filename, G4);
fprintf('\n');

% H is the header, S1 is the sequence of file 1, S2 is the sequence of file 2

H  = '[---[[[((((<FFFFFFFFF>))))][(((((-((((([[{I-IIIIII(([{I-I(-(({II(((--[[(((--((({I(--((((((((((-(<FF>))))))))-))))II}))))))----][((((---[(--[((((((((--(((--(((((((([[(-((((((((((((((<FF>)))))))))))-)))-)][((---((((<FF>))))--))]][((((((----(((<FF>)))))))))]))))))))-)))))))))))][((({I-III((--((((((((<FFFFFFF>))))))))))II--III})))])][(((((((((<FF>)))))---))))]))))]][[((((((((<FFF>)-)--))))))][((((((<FF>))))))]]--)))I})))II--I}][(-(({I--I((((((<FF>)))))-)I})))]-))IIIIIII}][([((((({I--I[(((((<FF>)))))]III})))))][{II-I[(((((((-(-(-{I((((((-(<FFFFFFF>)))-))))II})--)--)))))))-]III}]-)]][(({III----I[((((((---(((<FFFF>)))---)--)))))]III}))]--)))))))))-)]---][[(---{III------II(-[[((((((((((-[(((((((-(({I(((((((-----((((((((<FFF>))))))))---)))))))III})))))))))][{I-I(((((((((--((((((([((-((((((((((((((-(<FFFF>)-))---))))))))))))-))------][((<FFFF>))]-)))))))-)))))-))))III}]))))))))))][(-(((((((---(((((((((<FFFFFFFFF>))))))))))))))))-)]-----][{IIIIIII([(((((((((((((<FFFFF>)))))))))))))][(((<FFFF>)))]--)IIIIIII}])III})][((-(((-{II(---(((((<FF>))))))I})))))]]-------][[[{II[(((-(((((((-{I(--([(((((([(((((((((((-{I-I[-[(((<FFFF>)----))------------][((((((-[{II---II((((((----[(((((((<FFFF>)))))))-][(((<FFFF>)))]-)))))-)II}][{I((-(((--((({I((((((---[(((((([(((<FFFF>)))-][(((((<FFFF>)))))]--))))))----][((((-[(((((((--(((--{III<FFFF>IIII})))---)))))))][{III-I((((((<FFF>)-)))))III--II}]-))))]---)--)))-))I--I}))))))))I-I}]-))))))]]I}---))--)))))))))][([(--(((((((((--((((--((((((((<FFFF>))))))))----))))--)))))))))----)][(((((((((((((<FFF>)---))))))))))))]---)]-))))))][({II((((((((((<FFFF>)--)))))))))I})])--)I---I}))))))))))]II-I}][(--(-(((-(--(((((((((-(((((((((((((((-((((({II((((--(((<FF>)))--))))I})))))))))))))))))))))))))))))---)-))-)-)--)]-][(((((((((<FFFF>)-))))))))]]------]';

S1 = 'UGGA---GAGU<UUGAUCCUG>GCUC--AGGGUGAACGC--{UGGCGGCGUG-{CCUAAGA{CAUGCA---AGU--CGU{GCG-GGCCGCGGGGUU<UU>ACUCCGUGGUCAGCG}GCGGAC----G-GGUGAGU-AAC-GCGUGGGUGACCUACCCGGAAGA--GGGGGACAACCCGGGG<AA>ACUCGGGCUAAUCCCCC--AUGU-GGAC<CC>GCCCCUUG---GGGUGUGUCCAAA<GG>GCUUUGCCC-GCUUCCGGAUGGGCCCGCGU--CCC{AUCAGCUAGUUGGUGGG<GUAAUGG>CCCACCAAGGCGACGAC}GGG-UA-GCCGGUCUG<AG>AGGAUGGCCGGC-CACA----GGGGCACU<GAG>ACACGGGCCCC--ACUCCU<AC>GGGAGG--CAGCAG}UUAGGAAU}C-UUCC{GCAAUGGGCG<CA>AGCCUGAC}GGA-GCGACGCCGC}U-U-GGAGG{AAGAAGCCCU<UC>GGGGUGUAA}ACUCC--{UG-AACCCGGG--A-C-{GAAACCC-C<CGACGAG>GGG-ACUGAC}G-GU---ACCGGGGUAAU}--A---GC{GCCGGCCAACUCCGUGCCAGC<AGCC>GCGGUAAUACGGAG-GGC}GC-GAGCGUUACCC-G--GAU--UCA-{CUGGGCGUAAAGG--GCGUGUAGGCG-GCCUGGGGCG{UCCCAUGUGAAAGACCACGGC<UCA>ACCGUGGGGGAGCGUGGGAUA}CGCUCAGGC--{UAGACGGUGGGAGAGGGUGGU-GGAAUUCCCGGAGUAGCGG<UGAA>AUGCGCAGAUACCGGGAGGAACGCCGAUG-GC<GAAG>GC-AGCCACCUGGUCCACCCGUGAC}-GCUGAGGCGCG-AAAGCGUGGGGAGCAAACCGG<AUUAGAUAC>CCGGGUAGUCCACGCC-C---UAAA-{CGAUGCGC-GCUAGGUCUCUGG<GUCUC>CUGGGGGCCGAAGC-UAA<CGCG>UUA-AGCGCGCCGC}-CUGG}G--GAGUACG{GCCGCAAGGCU<GA>AACUCAA}AGGAA---UUGACGG---{GG-GCCCGCACAAGC{GGUGG-AGCAUG-UGGUUUAAUUC-{GAAGC-AAC<GCGA>AGAACCUUACCAGGCCUUGA-CAUGCU--{AGGGAAC-CCGGGUGAA-AGCCUGG<GGUG>CCCCGCG-A--GG<GGAG>CC---CUAGC--AC}--{AGGUGCUGCAUG{GCCGUCGUCA-GCUCGU-GCC<GUGA>GGUGU-UGGGU<UAAG>UCCCG-CAACGAGCGCAAC-CCCCG-CCGUUAGUUGCCAG{CGG<UUCG>GCCG}GGCACUCUAACGG--{GACUGCCCGC-<GAA>-AGCGGGAGGAAGG}--AGGG-GACGACGUCUGGUCAG}CAUGGCCCUUA}-CGGCCUG--G}-GCGACACACGUGCUA--C-AAUGCCCACUACAAAGCGAUGCCACCCG<GCAA>CGGGGAGCUAAUCGCAAAAAGGUGGGCCCAGUU-CGGAUUGGGGUCU<GCA>ACCCGACCCCAUGAAG-CCGG-AAUCGCU--A{GUAAUCGCGGAU<CAGC>CAUGCCGCGGUGA}A-U-ACGUUCC}CGGGCCUUGUACACA}--CCGCCCGUCAC-GCCAUGGGAGCGGGCUCUACCC-GAAGUCGC{CGGGAGCCUAC<GG>GCAGGCGCCG}AGGGUAG-GGCCCGUGACUGGGGCGAAGU-CGUAACAAGGUAG-CU-GUACCGGAA<GGUG>CGGCUGGAUC--ACUUUCU';

S2 = 'UGAA---GAGU<UUGAUCAUG>GCUC--AGAUUGAACGC--{UGGCGGCAGG-{CCUAACA{CAUGCAA--GUCGAACG{GUAACAGGAAGAAG-C<UU>GCUUCUUUGCUGACG}AGUGGCGGACG-GGUGAGU-AAU-GUCUGGGAA-ACUGCCUGAUGGA--GGGGGAUAACUACUGG<AA>ACGGUAGCUAAUACCGC--AUAACGUCG<CA>AGAC-CAA-A-GAGGGGGA--CCU<UC>GGGCCUCUU-GCCAUCGGAUGUGCCCAGAU--GGG{AUUAGCUAGUAGGUGGG<GUAACGG>CUCACCUAGGCGACGAU}CCC-UA-GCUGGUCUG<AG>AGGAUGACCAGC-CACA----CUGGAACU<GAG>ACACGGUCCAG--ACUCCU<AC>GGGAGG--CAGCAG}UGGGGAAU}A-UUGC{ACAAUGGGCG<CA>AGCCUGAU}GCA-GCCAUGCCGC}G-U-GUAUG{AAGAAGGCCU<UC>GGGUUGUAA}AGUAC--{UUUCAGCGGGGAGGAAG{GGAGUAAAG<UUAAUAC>CUUUGCUCAU}UGACGUUACCCGCAGAAG}-AA---GC{ACCGGCUAACUCCGUGCCAGC<AGCC>GCGGUAAUACGGAG-GGU}GC-AAGCGUUAAUCGG-AAUU--ACUG{GGCGUA---AAGC--GCACGCAGGCG-GUUUGUUAAG{UCAGAUGUGAAAUCCCCGGGC<UCA>ACCUGGGAACUGCAUCUGAUA}CUGGCAAGC--{UUGAGUCUCGUAGAGGGGGGU-AGAAUUCCAGGUGUAGCGG<UGAA>AUGCGUAGAGAUCUGGAGGAAUACCGGUG-GC<GAAG>GC-GGCCCCCUGGACGAAGACUGAC}-GCUCAGGUGCG-AAAGCGUGGGGAGCAAACAGG<AUUAGAUAC>CCUGGUAGUCCACGCCGU-AAACGA-{UGUCGACU-UGGAGGUUGUGCC<CUUGA>GGCGUGGCUUCCGG-AGC<UAAC>GCG-UUAAGUCGAC}-CGCC}U--GG-GGA-{GUAC--GGCCG<CA>AGGUUAA}AACUC--AAAUGAAU---{UG-ACGGGGGCCCGC{ACA-A-GCGGUG-GAGCAUGUGGUU{UAAUU-CGA<UGCA>A--CGCG--AAGAACCUUAC-CUGGUCU-{UGAC-AUCCACGGAA---GUUUUCA<GAGA>UGAGAAUGU-GCC<UUCG>GGA-ACCGUGAGAC}--{AGGUGCUGCAUG{GCUGUCGUCA-GCUCGU-GUU<GUGA>AAUGU-UGGGU<UAAG>UCCCG-CAACGAGCGCAAC-CCUUA-UCCUUUGUUGCCAG{CGG<UCCG>GCCG}GGAACUCAAAGGA--{GACUGCCAGUG<AUA>A-ACUGGAGGAAGG}-UGGGG--AUGACGUCAAGUCAU}CAUGGCCCUUA}-CGACCAG--G}GCUAC-ACACGUGCUA--C-AAUGGCGCAUACAAAGAGAAGCGACCUC<GCGA>GAGCAAGCGGACCUCAUAAAGUGCGUCGUAGUC-CGGAUUGGAGUCU<GCA>ACUCGACUCCAUGAAG--UCG-GAAUCGCU-A{GUAAUCGUGGAU<CAGA>A-UGCCACGGUGA}A-UACGU-UCC}CGGGCCUUGUACACA}--CC-GCCCGUCACACCAUGGGA-GUGGGUUGCAAAAGAAGUAGG{UAGCUUAACCU<UC>GGG-AGGGCG}CUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAG-GU-AACCGUAGG<GGAA>C-CUGCGGUU---GGAUCA';

%SeqAlignment = pInferAlignmentFromSequenceAlignment(File,S1,S2);

[cti1,itc1] = pColumnToIndex(S1);
[cti2,itc2] = pColumnToIndex(S2);

AMatrix2 = sparse(zeros(length(itc1),length(itc2)));

for c = 1:length(S1),
  if ismember(S1(c),'ACGUacgu') && ismember(S2(c),'ACGUacgu'),
    AMatrix2(cti1(c),cti2(c)) = 1;
  end
end

fprintf('%s\n',H);
fprintf('%s\n',S1);
fprintf('%s\n',S2);

Text{1} = '';

for c = 1:length(S1),
 if ismember(H(c),'<>{}'),
  Text{1} = [Text{1} sprintf(' ')];
 elseif S1(c) == '-' && S2(c) == '-',
  Text{1} = [Text{1} sprintf(' ')];
 elseif S1(c) == '-',
  y = find(AMatrix1(:,i2));
  if isempty(y),
    Text{1} = [Text{1} sprintf('.')];                     % neither claims an alignment
  else
    d = abs(min(y-i1));
    if d < 10,
      Text{1} = [Text{1} sprintf('%d',d)];
    else
      Text{1} = [Text{1} sprintf('X')];
    end
  end
 elseif S2(c) == '-',
  x = find(AMatrix1(i1,:));
  if isempty(x),
    Text{1} = [Text{1} sprintf('.')];
  else
    d = abs(min(x-i2));
    if d < 10,
      Text{1} = [Text{1} sprintf('%d',d)];
    else
      Text{1} = [Text{1} sprintf('X')];
    end
  end
 else                                 % sequence alignment claims an alignment
  i1 = cti1(c);
  i2 = cti2(c);
  if AMatrix1(i1,i2) == 1,            % 3D alignment has these
    if AMatrix2(i1,i2) == 1,          % and sequence alignment does too
      Text{1} = [Text{1} sprintf('=')];
    else
      Text{1} = [Text{1} sprintf('?')];                   % this should never appear
    end
  else                                % 3D alignment does not have these
    x = find(AMatrix1(i1,:));
    y = find(AMatrix1(:,i2));
    if isempty(x) && isempty(y),      % neither aligned in 3D alignment
      Text{1} = [Text{1} sprintf('.')];                   % no complaint    
    else
      if isempty(x),
        x = 10000;                    % make x non-empty
      else
        x = x(1);                     % take just one element
      end
      if isempty(y),
        y = 10000;                    % make y non-empty
      else
        y = y(1);                     % take just one element
      end
      d = abs(min(min(x-i2),min(y-i1)));   % how many columns off are we?
      if d < 10,
        Text{1} = [Text{1} sprintf('%d',d)];
      else
        Text{1} = [Text{1} sprintf('X')];
      end
    end
  end
 end
end  


fprintf('%s\n', Text{1});
fprintf('\n');
fprintf('The line below the sequence alignment compares the sequence alignment to the 3D alignment\n');
fprintf('= means that both alignments agree at this position\n');
fprintf('. means that the 3D alignment makes no claim about this position\n');
fprintf('1 means that the sequence alignment is off by one position right or left from the alignment claimed in the 3D alignment, similarly for 2, 3, ...\n');
fprintf('X means that the difference is more than 9\n');
