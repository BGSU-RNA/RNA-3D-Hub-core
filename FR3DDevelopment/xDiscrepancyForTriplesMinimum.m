% xDiscrepancyForTriplesMinimum(Model,Cand) considers all possible matchings between two triples and returns the minimum discrepancy

% Model and Cand are two lists of indices or nucleotide numbers

function [Disc] = xDiscrepancy(File1,Model,File2,Cand,LocationWeight,AngleWeight)

%figure(1)
%clf
%zDisplayNT(File1,Model);
%figure(2)
%clf
%zDisplayNT(File2,Cand);

% if File1 is a text string (filename), load the file

if strcmp(class(File1),'char'),
  File1name = File1;
  File1 = zAddNTData(File1name,2);
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(Model),'char'),
  Model = {Model};
end

if strcmp(class(Model),'cell'),
  Model = File1.NT(zIndexLookup(File1,Model));
else
  Model = File1.NT(Model);
end

% if File2 is a text string (filename), load the file and display

if strcmp(class(File2),'char'),
  File2name = File2;
  File2 = zAddNTData(File2name,2);
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(Cand),'char'),
  Cand = {Cand};
end

if strcmp(class(Cand),'cell'),
  Cand = File2.NT(zIndexLookup(File2,Cand));
else
  Cand = File2.NT(Cand);
end

if length(Model) ~= length(Cand),           % motif sizes must be the same
  Disc = [];
  R    = eye(3);
else

L = length(Cand);

A = zeros(L,1);                             % rotation angles for bases

if nargin < 5,
  LocationWeight = ones(1,L);
else
  LocationWeight = length(LocationWeight) * LocationWeight / sum(LocationWeight);
end

if nargin < 6,
  AngleWeight    = ones(1,L);
end

% ------------------------------ Calculate discrepancy

F1.NT = Model;
F2.NT = Cand;

d(1) = xDiscrepancyForTriples(F1,1:3,F2,[1 2 3]);
d(2) = xDiscrepancyForTriples(F1,1:3,F2,[2 1 3]);
d(3) = xDiscrepancyForTriples(F1,1:3,F2,[1 3 2]);
d(4) = xDiscrepancyForTriples(F1,1:3,F2,[3 2 1]);
d(5) = xDiscrepancyForTriples(F1,1:3,F2,[3 1 2]);
d(6) = xDiscrepancyForTriples(F1,1:3,F2,[2 3 1]);

Disc = min(d);

end
