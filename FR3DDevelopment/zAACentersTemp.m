

for f = 1:length(File),

flag = 0;

File(f).Filename

for i = 1:length(File(f).AA),
  if length(File(f).AA(i).Center) < 3,
    flag = 1;
    File(f).AA(i).Center = mean(File(f).AA(i).Loc,1);
  end
end

if flag == 1,
  zSaveNTData(File(f));
end

end
