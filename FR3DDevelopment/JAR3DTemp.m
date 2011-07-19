
clear java

A = Alignment.JAR3D(pwd,'16S_sequence_from_2avy_1j5e.fasta','16S_from_2AVY.txt',2,0,10);

Al = Alignment.getAlignment(A)

Al.get(0)
Al.get(1)
Al.get(2)
Al.get(3)
Al.get(4)

