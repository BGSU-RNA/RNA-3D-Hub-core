
% File = zAddNTData('1s72');
 File = zAddNTData('2zjr');
% File = zAddNTData('3G78');

Verbose = 1;

for f = 1:length(File),

Beta = [];
for n = 1:length(File(f).NT),
  i = find(File(f).NT(n).Beta < Inf);
  B = File(f).NT(n).Beta(i);
  Beta = [Beta; B];
  BetaRange(n) = max(B) - min(B);
end

figure(1)
clf
hist(Beta,30);
title(['Histogram of all atom beta factors ' File(f).Filename]);

saveas(gcf,['Beta factors ' File(f).Filename '.png']);

if Verbose > 1,

figure(2)
clf
hist(BetaRange,30) 
title('Range of beta factors within one nucleotide') 

for n = 1:length(File(f).NT),
  i = find(File(f).NT(n).Beta(1:12) < Inf);
  B = File(f).NT(n).Beta(i);
  BetaRange(n) = max(B) - min(B);
end

figure(3)
clf
hist(BetaRange,30) 
title('Range of beta factors within one sugar') 

for n = 1:length(File(f).NT),
  i = find(File(f).NT(n).Beta(13:end) < Inf);
  i = i + 12;
  B = File(f).NT(n).Beta(i);
  BetaRange(n) = max(B) - min(B);
end

figure(4)
clf
hist(BetaRange,30) 
title('Range of beta factors within one base') 

end

pause

end
