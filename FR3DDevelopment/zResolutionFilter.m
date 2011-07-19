% zResolutionFilter returns only files with resolution at or below ResLimit

function [File] = zResolutionFilter(File,ResLimit)

Res = zeros(1,length(File));

for f = 1:length(File),
  if ~isempty(File(f).Info.Resolution),
    Res(f) = File(f).Info.Resolution;
  end
end

i = find(Res <= ResLimit);

File = File(i);
