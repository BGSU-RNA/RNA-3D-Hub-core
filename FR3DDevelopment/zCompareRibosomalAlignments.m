
% Thermus first, E. coli second!

if 0 > 1,
  File(1) = zAddNTData('1j5e');
  File(2) = zAddNTData('2qan');
end

zCompareAlignmentToStructures(File,'Blake Ecoli TTh 16S michigan rdp_download_2seqs.fasta','Michigan16S')

zCompareAlignmentToStructures(File,'Blake Ecoli TTh 16S greengenes_export_0.fasta','Greengenes16S')

zCompareAlignmentToStructures(File,'Blake Ecoli Tth 16S arb-silva.de_2010-10-29_id5294.fasta','Silva16S')

zCompareAlignmentToStructures(File,'JAR3D Ecoli Tth 16S 2010-10-29.fasta','JAR3D')

