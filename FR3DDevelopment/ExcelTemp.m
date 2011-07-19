Excel = actxserver('Excel.Application');
Excel.Workbooks.Open('C:\Documents and Settings\zirbel\My Documents\FR3D\temp.xls');


NumbCol  = {'C','G','K','O'};
PercCol  = {'D','H','L','P'};
LeftCol  = {'E','I','M','Q'};
RightCol = {'F','J','N','R'};

for cc = 1:12,
  d = CI{cc};                            % data values to fill in
  for i = 1:4,                           % run through all pairs
   for j = 1:4,               
    pc = i + (j-1)*4;                    % paircode
    
    if any(cc == [3 4 5 6 9 10 11 12]) || (i <= j),

      Range = Excel.Range([NumbCol{j} num2str(4+20*(cc-1)+4*(i-1))]);
      Range.Value = d(pc,1);

      Range = Excel.Range([PercCol{j} num2str(4+20*(cc-1)+4*(i-1))]);
      Range.Value = d(pc,2);

      Range = Excel.Range([LeftCol{j} num2str(4+20*(cc-1)+4*(i-1))]);
      Range.Value = d(pc,3);

      Range = Excel.Range([RightCol{j} num2str(4+20*(cc-1)+4*(i-1))]);
      Range.Value = d(pc,4);

    end
   end
  end

  Range = Excel.Range(['S' num2str(4+20*(cc-1)+4*(i-1))]);
  Range.Value = length(BPCountsFromSequences{cc}(1,:));

end

Excel.Visible=1;
