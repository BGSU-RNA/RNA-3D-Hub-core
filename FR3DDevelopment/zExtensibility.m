


Extensible = zeros(1,length(File.NT));

E = fix(abs(File.Edge));
N = (E > 0) .* (E < 13) .* (File.Crossing > 1);     % non-nested basepairs only
E = E .* (E == 1) .* (File.Crossing == 0);    % nested cWW basepairs only
B = File.BasePhosphate;
B = (B < 100) .* (File.Crossing > 1);

% N = N + B;

for i = 1:length(File.NT),
  if Extensible(i) == 0,
    j = find( E(i,:) > 0 );
    if length(j) > 0,
      j = max(j);
      if j > i,
        if sum(sum(N(i:j, [1:(i-1) (j+1):length(File.NT)]))) == 0,
          Extensible(1,i) = 1;
          Extensible(1,j) = 1;
          Extensible(1,i:j) = ones(1,length(i:j));
        end
      end
    end
  end
end

i = find(Extensible == 1);

for ii = 1:length(i),
  fprintf('Extensible: %s%4s\n', File.NT(i(ii)).Base, File.NT(i(ii)).Number);
end