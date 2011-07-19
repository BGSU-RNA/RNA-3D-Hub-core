% zExemplarFrequencyCalculation retrieves every exemplar combination and calculates the IDI between it and every other pair in the family, saving the IDIs in the variable ExemplarIDI in PairExemplars.mat

Letters = 'ACGU';

if ~exist('Exemplar'),
  load PairExemplars
end

clear ExemplarFreq

for Class = [1:15],
fprintf('Calculating basepair frequencies for %s\n', zEdgeText(Class));
  for Code1 = 1:4,
    for Code2 = 1:4,

      [NT1,NT2,E] = zGetExemplar(Class,Code1,Code2);
      if ~isempty(E.Filename),
        Count(Code1,Code2) = E.Count;
      else
        Count(Code1,Code2) = 0;
      end
    end
  end

Count

  ExemplarFreq{Class} = Count;

end

save(['FR3DSource' filesep 'PairExemplars'],'Exemplar','ExemplarIDI','ExemplarFreq'); % Matlab version 7 only
save PairExemplars_Version_6.mat Exemplar ExemplarIDI ExemplarFreq -V6 % for compatibility with older versions
