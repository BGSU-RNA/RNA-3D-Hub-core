
Letter = 'ACGU';

BP1 = 'cWW';
BP2 = 'cHW';

    for a = 1:4,
      for b = 1:4,
        for c = 1:4,
figure(1)
clf
          F = zMakeTriple(BP1,BP2,a,b,c);
          Filename = ['Triple_' BP1 '_' BP2 '_' Letter(a) Letter(b) Letter(c) '.pdb'];
          if ~isempty(F),
            title(strrep(Filename,'_','-'));
            drawnow
            pause
            zWritePDB(F,Filename);
          end
        end
      end
    end

