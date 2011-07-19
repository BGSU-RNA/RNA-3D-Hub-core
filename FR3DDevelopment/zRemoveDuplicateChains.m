% Some PDB files contain multiple copies of the same molecule, not
% interacting with anything else.  Others, using Biological Unit
% Coordinates, have a chain interacting with another copy of itself.

function [File] = zRemoveDuplicateChains(File,Verbose)

