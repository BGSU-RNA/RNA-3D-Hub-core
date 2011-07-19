% zAttachAlignment reads the spreadsheet StructureToAlignmentMap to see if there is an entry for File, and if so, stores the variable FastaCol in File.  FASTAFiles is a cell array indicating the .fasta file for each structure.  It is optional.  LineNumber indicates which line of each fasta file is to be used to make the correspondence to the 3D structure.

% A current problem is that there may be multiple entries for the PDB ID in the structure to alignment Excel file.  Which one to use?

% File = zAddNTData('3i8h');
% F = zAttachAlignment(File,2,{'16S_Bacterial_Stombaugh_et_al_Sup_Mat_S2.fasta'},{'A'},0)
% F = zAttachAlignment(File,2,{'16S_Bacterial_Stombaugh_et_al_Sup_Mat_S2.fasta'},{'A'},312)

function [File] = zAttachAlignment(File,Verbose,FASTAFile,FASTAChain,LineNumber)

if nargin < 2,
  Verbose = 1;
end

if nargin < 3,
  % [n,t,r] = xlsread(['Alignments' filesep 'StructureToAlignmentMap.xls']);
  [n,t,r] = xlsread([DropboxRoot filesep 'Alignments\StructureToAlignmentMap.xls']);  % search current path

  for f = 1:length(File),
    p = find(ismember(upper(t(:,1)),upper(File(f).Filename)));  % find line(s)
    if length(p) > 0,                                 % if one or more,
      for a = 1:length(p),                          % 23S and 5S, difft aligs
          FASTA = zReadFASTA(r{p(a),3});
          if Verbose > 0,
            fprintf('Read alignment %s for %s\n', r{p(a),3}, File(f).Filename);
          end

          Chain = r{p(a),2};
          Entry = r{p(a),4};

          if r{p(a),4} > 0,
            File(f) = zAlignToFASTA(File(f),r{p(a),2},FASTA,r{p(a),4},Verbose);
          else
            File(f) = zAlignToFASTA(File(f),r{p(a),2},FASTA,0,Verbose);
          end

          if Verbose > 0,
            fprintf('Incorporated sequence data with %s\n', File(f).Filename);
          end
      end
    else
      if Verbose > 0 && ~isempty(File(f).NT),
        fprintf('zAttachAlignment: No alignment found for %s.\n', File(f).Filename);
      end
    end
  end

else

  

  for f = 1:length(File),
    FASTA = zReadFASTA(FASTAFile{f});
    File(f) = zAlignToFASTA(File(f),FASTAChain{f},FASTA,LineNumber(f),Verbose);
  end

end
