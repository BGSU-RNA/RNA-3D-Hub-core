
% pInferVsActual3DStructure aligns the inferred 3D structure to the 3D
% structure inferred from an alignment to a structure-based model

% JH and JS from original 1j5e model:
% JAR3D header:

JH = '[---[[[(((((<FFFFFF>)-))))][(((((-((((([[[{I-IIIIII(([{I-I[-(({II(((-[[((((----(((((((((((((<FFFFF>)))))-)))))))))-)))------][((((---[(--[((((((((-(((--(((((((([[((((((((((((((((<FF>)))))))))))--)))))][(({III((((((((-(<FFFF>)))))))))I}))]][((((<FFF>))))]))))))))-)))))))))))][((({I-I{I(((--((((((((<FFFFFFF>)))))))))))I--I}II})))])][(((((((((<FF>)))))---))))]))))]][[((((((((<FFF>)-)--))))))][((((((<FF>))))))]]--)))I}))]II--I}][((-(({II-II[((((<FFFF>))))]II}))))]))IIIIIII}-][((((({I-III(((((<FF>)))))IIII})))))]][{II((-(((((((---(---((--(((((<FFF>)))))--))----)----)))))))))II}]][(--(({III---II[((((((---(((<FFFF>)))---)--)))))]III}))-)])))))))))-)]---][[(---{III---II(-[[((((((((((-[(((((((-(({I((((((((----(((((((<FFFFF>)))))))--))))))))III})))))))))][{I-I(((((((((--((((((([((((((((((((({III(((<FFFF>)))II--II})))))))))))))---][(--((<FFFF>)))])))))))-))))))-)))III}]))))))))))][(-(((((((---(((((((((<FFFFFFFFF>))))))))))))))))---)]-----------][{IIIIIII([((((((((((<FFFFF>))))))))))][(((<FFFF>)))]--)IIIIIII}])III})][((((({II(-(((((<FF>))))))I})))))]]-------][[[{II[(((-(((((((-{I((([(((((([(((((((((((------[((<FFFFFFFF>))---------][((((((({I[[{II--II((((((----[((((((<FF>))))))--][((((<FFFF>))))])))))-)II}][(((-(((--{III((((((----[(((((([(((<FFFF>)))-][((((((<F>)-)))))]--))))))---][((((--[(((((-(-(-(---((((<FF>)-))))----))-)))))][((--(((((((<FF>-)))))))-))]----))))]-----)))-)))I--III}))))))]]I}-)))))))]---))-)))))))))][([(--(((((((((--((((--(({II((((<FFFF>))))I}-))----))))--)))))))))----)][(((((((((((((<FFF>)---))))))))))))]--)]-))))))][({II((((((((((<FFFFF>))))))))))I})])-))I--I}))))))))))]II-I}][(-(((((((-((((((((((((((((((((((((-((((((((------((((<FF>)-))))-)))))))))))))))))))))))))))))))--))))-)))--)]-][(((((((((<FFFF>)))))))))]]-----]';

% JAR3D sequence:

JS = 'UGAA---GAGUU<UGAUCA>UGGCUC--AGAUUGAACGC---{UGGCGGCAGG-{CCUAACA{CAUGCA--AGUCGAACGGUAACAGGAAGA<AGCUU>GCUUCUUUGCUGACGAGUGGCGGACG-GGUGAGU-AAU-GUCUGGGAAACUGCCUGAUGGA--GGGGGAUAACUACUGG<AA>ACGGUAGCUAAUACCGCA---U{AACGUCGCAAGAC<CAAA>GAGGGGGACC}U----UCGG<GCC>UCUU-GCCAUCGGAUGUGCCCAGAU--GGG{AUU{AGCUAGUAGGUGGG<GUAACGG>CUCACCUAGGCGACG}AU}CCC-UA-GCUGGUCUG<AG>AGGAUGACCAGC-CACA----CUGGAACU<GAG>ACACGGUCCAG--ACUCCU<AC>GGGAGG--CAGCAG}UGGGGAAU}--AUUGC{ACAAU-GGGC<GCAA>GCCUGAU}GCAG-CCAUGCCGC}GU-GUAUG{AAGAAGGCCU<UC>GGGUUGUAA}AGUAC---{UUUCAGCGGGGAGGAAGGGAGUAAAGUU<AAU>ACCUUUGCUCAUUGACGUUACCCGCAGAA}---GAAGC{ACCGGCUAACUCCGUGCCAGC<AGCC>GCGGUAAUACGGAG-GGU}GCAA-GCGUUAAUCGG-AAUU--ACUG{GGCGUAAAGC--GCACGCAGGCG-GUUUGUUAAG{UCAGAUGUGAAAUCCCCGGG<CUCAA>CCUGGGAACUGCAUCUGAUA}CUGGCAAGC--{UUGAGUCUCGUAGAGGGGGGU-AGAAUUCCAGGUG{UAGCGG<UGAA>AUGCGUAGA}GAUCUGGAGGAAUACCG-GUGGC<GAAG>GCG-GCCCCCUGGACGAAGACUGAC}-GCUCAGGUGCG-AAAGCGUGGGGAGCAAACAGG<AUUAGAUAC>CCUGGUAGUCCACGCCGUAA-ACGAUGUCGACU-{UGGAGGUU-GUGCCCUUG-<AGGCG>-UGGCUUCCGG-AGC<UAAC>GCG-UUAAGUCGAC}-CGCC}U--GGGGA{GUACGGCCG<CA>AGGUUAA}AACUC--AAAUGAAU---{UG-ACGGGGGCCCGC{ACAA-GCGGUG-GAGCAUGUGGUUUAAUU-CG<AUGCAACG>CGAAGAACCUUA-CCUGGUC{U--{UGACAUCCACGGAAGU-UUUCAG<AG>AUGAGAAUG-UGCC<UUCG>GGAA-CCGUGAGAC}--AGGUGCUGC{AUGGCUGUCGUCA-GCUCGU-GUU<GUGA>AAUGU-UGGGUU<A>AGUCCCG-CAACGAGCGCAA-CCCUUA-UCCUUUGUUGCCAGCGGU<CC>GGCCGGGAACUCAAAGGA--GACUGCCAGUG<AU>AAACUGGAGGA-AGGUGGGG-AUGACGUCAAGUCAUCAU}GGCCCU-UA}CGACCAGG-GCUACACACGUGCUA--C-AAUGGCGCAUACAAAGAGAAGC{GACCUC<GCGA>GAGCA}AGCGGACCUCAUAAAGUGCGUCGUAGUC-CGGAUUGGAGUCU<GCA>ACUCGACUCCAUGAAG-UCG-GAAUCGCU-A{GUAAUCGUGGAU<CAGAA>UGCCACGGUGA}A-UACGUUCC}CGGGCCUUGUACACA}--CCGCCCGUCACACCAUGGGAGUGGGUUGCAAAAGAAGUAGGUAGCUUAACCUU<CG>GGAGGGCGCUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAG-GU-AACCGUAGG<GGAA>CCUGCGGUU--GGAUCA'; %  16S from 2avy -1797.3639874395858 



EH = '[---[[[((((<FFFFFFFFF>))))][(((((-((((([[{I-IIIIII(([{I-I(-(({II(((--[[(((--((({I(--(((((((((((<FF>))))))))-))))II}))))))----][((((---[(--[((((((((-(((--(((((((([[(-((((((((((((((<FF>)))))))))))-)))-)][((---((((<FF>))))-))]][((((((--(((<FF>)))))))))]))))))))-)))))))))))][((({I-III((--((((((((<FFFFFFF>))))))))))II--III})))])][(((((((((<FF>)))))---))))]))))]][[((((((((<FFF>)-)--))))))][((((((<FF>))))))]]--)))I})))II--I}][((-(({I--I((((((<FF>)))))-)I}))))]))IIIIIII}][([((((({I--I[(((((<FF>)))))]III})))))][{II-I[(((((((((-(-{I((((((-(<FFFFFFF>)))-))))II})--))-)))))))-]III}]-)]][(({III----I[((((((---(((<FFFF>)))---)--)))))]III}))]--)))))))))-)]---][[(---{III---II(-[[((((((((((-[(((((((-(({I(((((((-----((((((((<FFF>))))))))---)))))))III})))))))))][{I-I(((((((((--((((((([((-((((((((((((((-(<FFFF>)-))---))))))))))))-))---][(--((<FFFF>)))])))))))-)))))-))))III}]))))))))))][(-(((((((---(((((((((<FFFFFFFFF>))))))))))))))))-)]-----][{IIIIIII([(((((((((((((<FFFFF>)))))))))))))][(((<FFFF>)))]--)IIIIIII}])III})][((((({II(-(((((<FF>))))))I})))))]]-------][[[{II[(((-(((((((-{I(-([(((((([(((((((((((-{I-I[-[(((<FFFF>)--))----------][((((((-[{II--II((((((--[(((((((<FFFF>)))))))-][(((<FFFF>)))]-)))))-)II}][{I((-(((--((({I((((((---[(((((([(((<FFFF>)))-][(((((<FFFF>)))))]--))))))----][((((-[(((((((--(((--{III<FFFF>IIII})))---)))))))][{III-I((((((<FFF>))))))I-II-II}]-))))]--)--)))-))I--I}))))))))I-I}]-))))))]]I}---))-)))))))))][([(--(((((((((--((((--((((((((<FFFF>))))))))----))))--)))))))))----)][(((((((((((((<FFF>)---))))))))))))]--)]-))))))][({II((((((((((<FFFF>)-)))))))))I})])--)I--I}))))))))))]II-I}][(-(-(((-(--((((((((((((((((((((((((-((((({II((((--(((<FF>)))-))))I})))))))))))))))))))))))))))))---)-))-)-)--)]-][(((((((((<FFFF>)))))))))]]-----]';

ES = 'UGAA---GAGU<UUGAUCAUG>GCUC--AGAUUGAACGC--{UGGCGGCAGG-{CCUAACA{CAUGCAA--GUCGAACG{GUAACAGGAAGAAGC<UU>GCUUCUUUGCUGACG}AGUGGCGGACG-GGUGAGU-AAU-GUCUGGGAAACUGCCUGAUGGA--GGGGGAUAACUACUGG<AA>ACGGUAGCUAAUACCGC--AUAACGUCG<CA>AGACCAA-A-GAGGGGGACCU<UC>GGGCCUCUU-GCCAUCGGAUGUGCCCAGAU--GGG{AUUAGCUAGUAGGUGGG<GUAACGG>CUCACCUAGGCGACGAU}CCC-UA-GCUGGUCUG<AG>AGGAUGACCAGC-CACA----CUGGAACU<GAG>ACACGGUCCAG--ACUCCU<AC>GGGAGG--CAGCAG}UGGGGAAU}--AUUGC{ACAAUGGGCG<CA>AGCCUGAU}GCAG-CCAUGCCGC}G-U-GUAUG{AAGAAGGCCU<UC>GGGUUGUAA}AGUAC--{UUUCAGCGGGGAGGAAG{GGAGUAAAG<UUAAUAC>CUUUGCUCAU}UGACGUUACCCGCAGAAG}-AA---GC{ACCGGCUAACUCCGUGCCAGC<AGCC>GCGGUAAUACGGAG-GGU}GC-AAGCGUUAAUCGG-AAUU--ACUG{GGCGUAAAGC--GCACGCAGGCG-GUUUGUUAAG{UCAGAUGUGAAAUCCCCGGGC<UCA>ACCUGGGAACUGCAUCUGAUA}CUGGCAAGC--{UUGAGUCUCGUAGAGGGGGGU-AGAAUUCCAGGUGUAGCGG<UGAA>AUGCGUAGAGAUCUGGAGGAAUACCG-GUGGC<GAAG>GCG-GCCCCCUGGACGAAGACUGAC}-GCUCAGGUGCG-AAAGCGUGGGGAGCAAACAGG<AUUAGAUAC>CCUGGUAGUCCACGCCGU-AAACGA-{UGUCGACU-UGGAGGUUGUGCC<CUUGA>GGCGUGGCUUCCGG-AGC<UAAC>GCG-UUAAGUCGAC}-CGCC}U--GGGGA{GUACGGCCG<CA>AGGUUAA}AACUC--AAAUGAAU---{UG-ACGGGGGCCCGC{ACAA-GCGGUG-GAGCAUGUGGUU{UAAUU-CGA<UGCA>ACGCGAAGAACCUUAC-CUGGUCU-{UGACAUCCACGGAA-GUUUUCA<GAGA>UGAGAAUGU-GCC<UUCG>GGA-ACCGUGAGAC}--{AGGUGCUGCAUG{GCUGUCGUCA-GCUCGU-GUU<GUGA>AAUGU-UGGGU<UAAG>UCCCG-CAACGAGCGCAAC-CCUUA-UCCUUUGUUGCCAG{CGG<UCCG>GCCG}GGAACUCAAAGGA--{GACUGCCAGUG<AUA>AACUGGAGGAAGG}-UGGGG-AUGACGUCAAGUCAU}CAUGGCCCUUA}-CGACCAG--G}GCUACACACGUGCUA--C-AAUGGCGCAUACAAAGAGAAGCGACCUC<GCGA>GAGCAAGCGGACCUCAUAAAGUGCGUCGUAGUC-CGGAUUGGAGUCU<GCA>ACUCGACUCCAUGAAG-UCG-GAAUCGCU-A{GUAAUCGUGGAU<CAGA>AUGCCACGGUGA}A-UACGUUCC}CGGGCCUUGUACACA}--CCGCCCGUCACACCAUGGGAGUGGGUUGCAAAAGAAGUAGG{UAGCUUAACCU<UC>GGGAGGGCG}CUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAG-GU-AACCGUAGG<GGAA>CCUGCGGUU--GGAUCA'; %  16S from 2avy

j = 1;            % position in Java sequence
e = 1;            % position in experimental sequence

[Jcti,Jitc] = pColumnToIndex(JS);
[Ecti,Eitc] = pColumnToIndex(ES);

Text{1} = '';
Text{2} = '';
Text{3} = '';
Text{4} = '';

for i = 1:length(Jitc),
  Text{1} = [Text{1} JH(Jitc(i))];
  Text{2} = [Text{2} JS(Jitc(i))];
  Text{3} = [Text{3} EH(Eitc(i))];
  Text{4} = [Text{4} ES(Eitc(i))];
end

fprintf('%s Header from Java parser\n', Text{1});
fprintf('%s Sequence from Java parser\n', Text{2});
fprintf('%s Sequence from 3D structure\n', Text{4});
fprintf('%s Header from 3D structure\n', Text{3});
