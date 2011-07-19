

% make a folder called Triples in the FR3D folder!
Letter = 'ACGU';

for i = 1:12,
  BP1 = zEdgeText(i);
  for j = [-3 -10 -9 -6 -5 -4 1:12],

    if j ~= 0 && abs(j) >= abs(i),

      BP2 = zEdgeText(j);
      if BP1(3) ~= BP2(2),

        fprintf('Making %s - %s triples\n', BP1, BP2);

        for a = 1:4,
          for b = 1:4,
            for c = 1:4,
              F = zMakeTriple(BP1,BP2,a,b,c);
              if ~isempty(F),

% display in a figure window

                figure(1)
                clf
                VP.Sugar = 1;
                VP.AtOrigin = 1;
                zDisplayNT(F,[2 1 3],VP)
                Title = '';
                Title = [Title F.NT(1).Base F.NT(1).Number '-'];
                Title = [Title F.NT(2).Base F.NT(2).Number '-'];
                Title = [Title F.NT(3).Base F.NT(3).Number ' '];
                Title = [Title strrep(BP1,' ','') '-' strrep(BP2,' ','')];
                title(Title);
                view(2)

                Filename = ['Triples' filesep 'Triple_' BP1 '_' BP2 '_' Letter(a) Letter(b) Letter(c) '.png'];

axis off

saveas(gcf,Filename,'png');

close all
                Filename = ['Triples' filesep 'Triple_' BP1 '_' BP2 '_' Letter(a) Letter(b) Letter(c) '.pdb'];


%cat(1,F.NT.Fit)
%mean(cat(1,F.NT.Fit))

%                zWritePDB(F,Filename,F.NT(1).Rot,mean(cat(1,F.NT.Fit)));
                zWritePDB(F,Filename,F.NT(2).Rot,F.NT(2).Fit(1,:));



%                pause
              end
            end
          end
        end
      end
    end
  end
end
