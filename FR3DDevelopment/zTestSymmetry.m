% zTestSymmetry(File) looks at cWW, tWW, cHH, tHH basepairs to see if there is a 3' versus 5' orientation effect

% File = zAddNTData('NonRedundant_2008_02_21_list');
% zCountInteractions(File,'Basepair_counts');

function [Count,CCount] = zTestSymmetry(File)

CI  = [];                                     % count interactions
CAS = [];                                    % count anti/syn

categories = [1 2 7 8];
Letter = 'ACGU';

for c = 1:length(categories),
 for a = 1:4,
  for b = 1:4,
   x = 0;
   y = 0;
   for f = 1:length(File),
    E = fix(File(f).Edge);

    [i,j,k] = find(E==categories(c));

    if length(i) > 0,
      Code1 = cat(1,File(f).NT(i).Code);
      Code2 = cat(1,File(f).NT(j).Code);

      m = find( (Code1 == a) .* (Code2 == b));
      x = x + sum(i(m) < j(m));
      y = y + sum(i(m) > j(m));

      xx = sum(i(m) < j(m));
      yy = sum(i(m) > j(m));

      % ------------- show particularly extreme structures

      if (a == 1) && (b == 3) && 2*binocdf(min(xx,yy),xx+yy,0.5) < 0.1,
        File(f).Filename
        File(f).Info
        cat(1,File(f).NT.Base)'
        cat(1,File(f).NT.Chain)'
        [xx yy]
      end

    end
   end

   if (x+y > 0) && (a~=b),
     fprintf('%s%s%s i<j %4d i>j %4d Total %5d  p-value %8.4f\n',Letter(a),Letter(b),zEdgeText(categories(c)),x,y,x+y,2*binocdf(min(x,y),x+y,0.5));
   end
  end
 end
end

