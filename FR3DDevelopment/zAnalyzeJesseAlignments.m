
[File1,i1,File2,i2] = zReadJesseAlignments('16S');

File(1) = File1;
File(2) = File2;

zAlignmentDiagram(File,i1',i2');

