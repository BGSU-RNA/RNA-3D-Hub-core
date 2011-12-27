% zRetrieveAlignment queries Blake Sweeney's database of alignments

% the text string used here are JSON arrays.  There are programs that can convert JSON arrays into Matlab cell array.  See:
% http://www.mathworks.com/matlabcentral/fileexchange/23393


T = '{ pdb:"2AW7", chain:"A", positions:[["10","11","12","13"],["20","21","23","24"]],sequence-ids:"all", source:"greengenes", version:"current"}';

T = '{ pdb:"2AW7", chain:"A", positions:["10","11","12","13"],sequence-ids:"all", source:"greengenes", version:"current"}';

s = urlread('http://rna.bgsu.edu/api/variations','POST',T);





T = '{ pdb:"2AW7", chain:"A", positions:[[["10","11","12","13"],["20","21","23","24"]]],sequence-ids:"all", source:"greengenes", version:"current"}';

s = urlread('http://rna.bgsu.edu/api/motifs','POST',T);

s = urlread('http://www-math.bgsu.edu/z');

