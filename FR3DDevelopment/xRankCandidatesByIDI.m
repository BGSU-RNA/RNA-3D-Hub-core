% xRankCandidatesByIDI(File,Model,Candidates) calculates the discrepancy for
% each candidate

function [IDI, Candidates] = xRankCandidatesByIDI(File,Cand,Verbose)

M1 = File(Cand(1,3)).NT(Cand(1,1));
M2 = File(Cand(1,3)).NT(Cand(1,2));

if nargin < 3,
  Verbose = 0;
end

if Verbose > 0,
  fprintf('Calculating IDI\n');
end

s = length(Cand(:,1));

IDI = zeros(65000,1);
Candidates  = uint16(zeros(65000,length(Cand(1,:))));

count = 0;

tic

if Verbose > 0,
  fprintf('Seconds remaining:');
end

for i=1:s,
  N1 = File(Cand(i,3)).NT(Cand(i,1));
  N2 = File(Cand(i,3)).NT(Cand(i,2));
  A = zIsoDiscrepancy(M1,M2,N1,N2);

  if A >= 0,
    count = count + 1;
    IDI(count,1) = A;
    Candidates(count,:)  = Cand(i,:);
  end

  if (mod(i,round(s/10)) == 0) && (Verbose > 0)
    fprintf(' %d', fix((s-i)*toc/i)); 
    drawnow
  end
end

if Verbose > 0,
  fprintf('\n');
end

Candidates  = Candidates(1:count,:);
IDI         = IDI(1:count,:);

[y,i]       = sort(IDI);                    % sort by discrepancy
Candidates  = Candidates(i,:);
IDI = IDI(i);

if Verbose > 1,
  fprintf('Calculating discrepancy took        %8.3f seconds\n',toc);
end

drawnow

