% dNeedlemanWunsch(seq1,seq2,p,d) aligns sequences seq1 and seq2 using probability p of base conservation and gap penalty d

function [matches,align1,align2,s1,s2] = zNeedlemanWunsch(seq1,seq2,gap,extendgap)

if nargin < 3,
  gap = 8;                              % Matlab default open gap penalty
end

if nargin < 4,
  extendgap = gap;
end

[Score, Alignment, Start] = nwalign(seq1,seq2,'Alphabet','NT','GapOpen',gap,'ExtendGap',extendgap);

i = find(Alignment(1,:) ~= '-');        % locations of letters in seq1
align1 = find(Alignment(2,i) == '|');        % aligned indices of seq1

i = find(Alignment(3,:) ~= '-');        % locations of letters in seq2
align2 = find(Alignment(2,i) == '|');        % aligned indices of seq2

s1 = Alignment(1,:);
s2 = Alignment(3,:);

matches = sum(seq1(align1) == seq2(align2));





return

seq1 = 'aacguuguggaa';
seq2 = 'accguaugugcaga';
[m,a1,a2,s1,s2] = zNeedlemanWunsch(seq1,seq2)

[Score, Alignment, Start] = swalign('uuuuaacguuguggaa','aacguaugugcaga','Alphabet','NT')

File = zAddNTData({'1s72','2aw4'});
seq1 = cat(2,File(1).NT.Base);
seq2 = cat(2,File(2).NT.Base);
[Score, Alignment, Start] = nwalign(seq1,seq2,'Alphabet','NT')
