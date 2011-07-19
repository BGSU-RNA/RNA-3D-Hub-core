% zLongestChain(File) returns the base sequence of the longest RNA chain in File along with the starting and ending indices

function [Seq,i1,i2] = zLongestChain(File)

if length(File.NT) > 0,

  Bases = cat(2,File.NT.Base);                % all bases
  Chain = cat(2,File.NT.Chain);               % all chain identifiers
  U = unique(Chain);                          % unique chain identifiers

  if length(U) > 1,                           % more than one chain

    maxlen = 0;                               % current longest length

    for u = 1:length(U),                      % loop through chains
      i = find(Chain == U(u));                % indices of this chain
      if length(i) > maxlen,
        Seq = Bases(i);                       % store bases of this chain
        i1 = i(1);                            % starting index of this chain
        i2 = i(end);                          % ending index of this chain
      end
    end

  else

    Seq = Bases;
    i1 = 1;
    i2 = length(File.NT);

  end

else

  Seq = '';
  i1 = 0;
  i2 = 0;

end


