% zSubFile extracts nucleotides specified in i

function [F] = zSubFile(File,i)

F = File;
F.NT = File.NT(i);
F.Edge = File.Edge(i,i);
F.Coplanar = File.Coplanar(i,i);
F.Covalent = File.Covalent(i,i);
F.Crossing = File.Crossing(i,i);
F.BasePhosphate = File.BasePhosphate(i,i);
F.Backbone = File.Backbone(i,i);
F.Range = File.Range(i,i);
F.NumNT = length(i);
