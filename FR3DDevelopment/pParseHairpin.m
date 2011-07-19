
function [P] = zParseHairpin(node,P,x,n,s,i,j)

% In the future, the probabilities should depend on x(i:j)!
% They should also depend on the hairpin type

if (j-i <= 4),
  switch node(n).subtype
   case 'GNRA'
     if (j-i == 3),
       if ((x(j) == 'A') | ((x(i) ~= 'G') & (x(j) == 'C'))) & ...
         (x(i+2) == 'A' | x(i+2) == 'G'),
           P(n,s).mp(i,j) = log(10/904);
       else 
           P(n,s).mp(i,j) = log(1/904);
       end
     end
   case '....'
     P(n,s).mp(i,j) = -1000;
     if (j-i == 3),
       if (x(i) == '.') & (x(i+1)=='.') & (x(i+2)=='.') & (x(j)=='.'),
         P(n,s).mp(i,j) = 1;
       end
     end
   otherwise
     P(n,s).mp(i,j)   = log(1/(4+16+64+256));
   end
end

P(n,s).sub(i,j)  = i;                     % start of subseq for next node
P(n,s).sub(j,i)  = j;                     % end of subsequence
P(n,s).next(i,j) = 0;                     % child state
