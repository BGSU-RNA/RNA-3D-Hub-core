% zCompareOrganismNames(in,jn) compares the first comma-separated names in in and jn, returning 1 if they are different and 0 if they are the same.  We assume that the longest RNA chain is listed first.

% We presume that the organism names are different, and look for evidence that they are the same
% Chains that are designated as "synthetic" or "engineered" will not be recognized as being the same as any other name.  They will not match "synthetic" or "Escherichia coli"
% Different strains of the same organism will (at least should) match

function [CompareOrg] = zCompareOrganismNames(in,jn)

in = strrep(in,'synthetic','');               % remove "synthetic"
jn = strrep(jn,'synthetic','');
in = strrep(in,'engineered','');               % remove "engineered"
jn = strrep(jn,'engineered','');

in = [',' in ','];                            % pad with commas
jn = [',' jn ','];

ic = strfind(in,',');                         % locations of commas
jc = strfind(jn,',');

i = 1;
j = 1;

CompareOrg = 0;                               % cannot say they are different

if ic(i+1)-ic(i) > 7 && jc(j+1)-jc(j) > 7,      % enough text to compare
  is = in((ic(i)+1) : (ic(i+1)-1));
  js = jn((jc(j)+1) : (jc(j+1)-1));
  [namematches,na,nb,nss,ntt] = dNeedlemanWunsch(is,js,0.99,2);
  namematches = sum(is(na) == js(nb));

  namelength = min(length(is), length(js));
  namelength = max(namelength,1);       % check bio source for ribosomes
                                        % some have similar sequences

  if namematches <= 0.5 * namelength,
    CompareOrg = 1;                     % definitely different
  end
end






if 0 > 1,

for i = 1:(length(ic)-1),
  for j = 1:(length(jc)-1),
    if ic(i+1)-ic(i) > 5 && jc(j+1)-jc(j) > 5,      % enough text
        is = in((ic(i)+1) : (ic(i+1)-1));
        js = jn((jc(j)+1) : (jc(j+1)-1));
        [namematches,na,nb,nss,ntt] = dNeedlemanWunsch(is,js,0.99,2);
        namematches = sum(is(na) == js(nb));

        namelength = min(length(is), length(js));
        namelength = max(namelength,1);       % check bio source for ribosomes
                                              % some have similar sequences

        if namematches > namelength / 2,
          CompareOrg = 0;                     % similar enough
        end
    end
  end
end

end
