% zParseFilenameNucleotideString(s) extracts the filename from the first 4 characters of s and the NTList from the remaining characters, after a space.  It could easily be improved to be more robust.

function [File,NTList,Indices] = zParseFilenameNucleotideString(s)

File = zAddNTData(s(1:4));

NTList = s(6:end);

Indices = zIndexLookup(File,NTList);
