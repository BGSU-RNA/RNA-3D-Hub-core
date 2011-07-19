% pClassify(S1,S2,N1,N2) displays the log probabilities for the same
% sequences parsed according to two different models.  Nodes to use are
% specified by N1 and N2.

function [void] = pClassify(S1,S2,N1,N2)

if nargin < 3,
  N1(1) = 1;
  N1(2) = length(S1(1).TraceInfo)+1;
  N2(1) = 1;
  N2(2) = length(S2(1).TraceInfo)+1;
end

clf
hold on

r = 1;

for k=1:length(S1),
  if N1(2) > length(S1(1).TraceInfo),
    P1 = S1(k).TraceInfo(1,N1(1)).mp+r*(rand-0.5);
  else
    P1 = S1(k).TraceInfo(1,N1(1)).mp-S1(k).TraceInfo(1,N1(2)).mp+r*(rand-0.5);
  end
  if N2(2) > length(S2(1).TraceInfo),
    P2 = S2(k).TraceInfo(1,N2(1)).mp+r*(rand-0.5);
  else
    P2 = S2(k).TraceInfo(1,N2(1)).mp-S2(k).TraceInfo(1,N2(2)).mp+r*(rand-0.5);
  end
  if S1(k).FastaNum(1) == 'A',
    plot(P1,P2,'r+');
  else
    plot(P1,P2,'k.');
  end
end

title('Archaeal versus bacterial model - Motif only','FontSize',14);
xlabel('Log probability according to Archaeal model','FontSize',14);
ylabel('Log probability according to Bacterial model','FontSize',14);


