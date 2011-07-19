
function [F] = zExtractChain(File,Chain)

i = zIndexLookup(File,Chain);
F = zSubFile(File,i);

