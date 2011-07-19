
function [void] = zStackedLogBar(bins,heights,color,orientation)

if nargin < 4,
  orientation = 'H';
end

if orientation == 'V',

for i = 1:length(bins)-1,
  for j = 1:(length(heights(1,:))-1),
    B = heights(i,j);
    T = heights(i,j+1);
    L = exp(bins(i));
    R = exp(bins(i+1));
    if color(j) == 'o',
      patch([L R R L L], [B B T T B], [255  	165  	0]/255);
    else
      patch([L R R L L], [B B T T B], color(j));
    end
    hold on
  end
end

else

for i = 1:length(bins)-1,
  for j = 1:(length(heights(1,:))-1),
    B = heights(i,j);
    T = heights(i,j+1);
    L = exp(bins(i));
    R = exp(bins(i+1));
    if color(j) == 'o',
      patch([B B T T B], [L R R L L], [255  	165  	0]/255);
    else
      patch([B B T T B], [L R R L L], color(j));
    end
    hold on
  end
end

end
