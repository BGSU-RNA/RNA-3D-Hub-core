
tic; [A1,A2,maximalClique,VMI,EM] = R3DAlign('1s72','10:530_0','3jyx','_4,1:420_5',0.7,5,40,'Greedy'); toc

tic; [B1,B2,maximalClique,VMI,EM] = R3DAlign('1s72','530:635_0','3jyx','420:650_5',0.7,5,70,'Greedy'); toc

tic; [C1,C2,maximalClique,VMI,EM] = R3DAlign('1s72','635:903_0','3jyx','650:940_5',0.7,3,70,'Greedy'); toc

tic; [D1,D2,maximalClique,VMI,EM] = R3DAlign('1s72','903:1300_0','3jyx','940:1376_5',0.7,5,40,'Greedy'); toc

tic; [E1,E2,maximalClique,VMI,EM] = R3DAlign('1s72','1300:1512_0','3jyx','1376:1613_5',0.7,3,40,'Greedy'); toc

tic; [F1,F2,maximalClique,VMI,EM] = R3DAlign('1s72','1512:1724_0','3jyx','1613:1881_5',0.7,4,50,'Greedy'); toc

tic; [G1,G2,maximalClique,VMI,EM] = R3DAlign('1s72','1724:2050_0','3jyx','1881:2349_5',0.7,3,60,'Greedy'); toc

tic; [H1,H2,maximalClique,VMI,EM] = R3DAlign('1s72','2050:2277_0','3jyx','2349:2510_5',0.7,3,80,'Greedy'); toc

tic; [I1,I2,maximalClique,VMI,EM] = R3DAlign('1s72','2277:2481_0','3jyx','2510:2815_5',0.7,3,80,'Greedy'); toc

tic; [J1,J2,maximalClique,VMI,EM] = R3DAlign('1s72','2660:2914_0','3jyx','2992:3394_5',0.7,3,80,'Greedy'); toc

AA = [A1; B1; C1; D1; E1; F1; G1; H1; I1; J1];
BB = [A2; B2; C2; D2; E2; F2; G2; H2; I2; J2];

tic; [AAA,BBB,maximalClique,VMI,EM] = R3DAlign('1s72','10:2914_0','3jyx','_4,1:3394_5',0.7,5,20,'Greedy',AA,BB); toc

tic; [AAA,BBB,maximalClique,VMI,EM] = R3DAlign('1s72','10:2914_0','3jyx','_4,1:3394_5',0.7,3,10,'Greedy',AA,BB); toc



% idea:  Align the first 200 nucleotides of each structure.
% find the last aligned pair, then align the *next* 200 of each after that
% in this way, one should be able to string together a set of local alignments
% some parts of it, especially at the ends of these arbitrary 200-nucleotide 
% segments, may not be aligned very well, but still, it would be more
% automated as a way of generating a seed alignment.  Then you can use
% the seed for the next round of alignments.
