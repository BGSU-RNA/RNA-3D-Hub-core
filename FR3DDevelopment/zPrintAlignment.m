
function [void] = zPrintAlignment(File,Alignment)

for c = 1:length(Alignment.Correspondence),
  Corr = Alignment.Correspondence(c);
  fprintf('%s', Alignment.Filename{1});
  for n = 1:length(Corr.File1),
    fprintf(' %s%s', File(1).NT(Corr.File1(n)).Base, File(1).NT(Corr.File1(n)).Number);
  end

  fprintf('  %s', Alignment.Filename{2});
  for n = 1:length(Corr.File2),
    fprintf(' %s%s', File(2).NT(Corr.File2(n)).Base, File(2).NT(Corr.File2(n)).Number);
  end

  fprintf('\n');
end


