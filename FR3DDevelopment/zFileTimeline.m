
% File = zAddNTData('NonRedundant_2008_02_21_list');

for f = 1:length(File),
  NumNT(f) = File(f).NumNT;
  if ~isempty(File(f).Info.ReleaseDate),
    Date(f) = datenum(File(f).Info.ReleaseDate, 'mm/dd/yyyy');
    fprintf('%10d %12s %4s %5d\n', Date(f), File(f).Info.ReleaseDate, File(f).Filename, File(f).NumNT);
  else
    Date(f) = Inf;    
    fprintf('No date available\n');
  end
end

[y,i] = sort(Date);

File = File(i);
Date = Date(i);
NumNT = NumNT(i);

Year = 1995 + (Date - datenum('01/01/1995','mm/dd/yyyy'))/365;

clf
subplot(2,1,1)
stairs(Year,1:length(File));

subplot(2,1,2)
stairs(Year,cumsum(NumNT));

fprintf('\n\n');

for f = 1:length(File),
  if ~isempty(File(f).Info.ReleaseDate),
%    fprintf('%10d %s\n', Date(f), File(f).Info.ReleaseDate);
  else
%    fprintf('No date available\n');
  end
end
