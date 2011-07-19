Filename = 'SARSA 16S.txt';

fid = fopen(Filename,'r');

if fid > 0

  L = 1;

  c = 1;

  while L > -1
    L = fgetl(fid);
    if L > -1                     % chokes on blank lines!
      if ~isempty(strfind(L,'Alignment of original')),
        L = -1;
      end
    end
  end

  First = [];
  Align = [];
  Second= [];

  L = 1;

  while L > -1,
    L = fgetl(fid);
    if isempty(strfind(L,'/pre')),
      First = [First  L];
      Align = [Align  fgetl(fid)];
      Second= [Second fgetl(fid)];
    else
      L = -1;
    end
  end

  fclose(fid);
end

First = First(22:(end-6));
Align = Align(22:end);
Second= Second(22:(end-6));

j = find(First ~= ' ');

First = First(j);
Align = Align(j);
Second = Second(j);

% [First; Align; Second]
