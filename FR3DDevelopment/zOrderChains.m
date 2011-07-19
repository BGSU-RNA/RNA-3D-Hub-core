% zOrderChains re-orders the nucleotides in File to put the longest chain first

function [File] = zOrderChains(File)

for f = 1:length(File),

  Chain = cat(2,File(f).NT.Chain);

  U = unique(Chain);

  Order = zeros(1,length(Chain));

  for u = 1:length(U),
    Order = Order + sum(Chain == U(u)) * (Chain == U(u));
  end

  Order(2,:) = 1:length(Chain);

  Order = Order';
  [M,i] = sortrows(Order, [-1 2]);

  if length(i) > 1,

    File(f).NT   = File(f).NT(i);
    File(f).Edge = File(f).Edge(i,i);
    File(f).Range = File(f).Range(i,i);
    File(f).BasePhosphate = File(f).BasePhosphate(i,i);
    if ~isempty(File(f).Distance),
      File(f).Distance = File(f).Distance(i,i);
    end
  end
end

