

load PairExemplars

  [s,t] = size(Exemplar);
  for a = 1:s,
    for b = 1:t,
      if ~isempty(Exemplar(a,b).Filename),
        FN = Exemplar(a,b).Filename;
        if ~isempty(strfind(FN,'*CURATED')) || ~isempty(strfind(FN,'*MODEL')),

clear NE

        NE.Filename = Exemplar(a,b).Filename;
        NE.Class = Exemplar(a,b).Class;

        if isfield(Exemplar(a,b).NT1,'Model'),
          NE.NT1.Model = Exemplar(a,b).NT1.Model;
        else
          NE.NT1.Model = 1;
        end
        NE.NT1.Base = Exemplar(a,b).NT1.Base;
        if isfield(Exemplar(a,b).NT1,'Unit'),
          NE.NT1.Unit = Exemplar(a,b).NT1.Unit;
        else
          NE.NT1.Unit = Exemplar(a,b).NT1.Base;
        end
        NE.NT1.Chain = Exemplar(a,b).NT1.Chain;
        NE.NT1.Number = Exemplar(a,b).NT1.Number;
        NE.NT1.Loc = Exemplar(a,b).NT1.Loc;
        NE.NT1.Sugar = Exemplar(a,b).NT1.Sugar;
        if isfield(Exemplar(a,b).NT1,'Beta'),
          NE.NT1.Beta = Exemplar(a,b).NT1.Beta;
        else
          NE.NT1.Beta = [];
        end
        NE.NT1.Center = Exemplar(a,b).NT1.Center;
        NE.NT1.Code = Exemplar(a,b).NT1.Code;
        NE.NT1.Rot = Exemplar(a,b).NT1.Rot;
        NE.NT1.Fit = Exemplar(a,b).NT1.Fit;
        NE.NT1.Syn = Exemplar(a,b).NT1.Syn;

        if isfield(Exemplar(a,b).NT2,'Model'),
          NE.NT2.Model = Exemplar(a,b).NT2.Model;
        else
          NE.NT2.Model = 1;
        end
        NE.NT2.Base = Exemplar(a,b).NT2.Base;
        if isfield(Exemplar(a,b).NT2,'Unit'),
          NE.NT2.Unit = Exemplar(a,b).NT2.Unit;
        else
          NE.NT2.Unit = Exemplar(a,b).NT2.Base;
        end
        NE.NT2.Chain = Exemplar(a,b).NT2.Chain;
        NE.NT2.Number = Exemplar(a,b).NT2.Number;
        NE.NT2.Loc = Exemplar(a,b).NT2.Loc;
        NE.NT2.Sugar = Exemplar(a,b).NT2.Sugar;
        if isfield(Exemplar(a,b).NT2,'Beta'),
          NE.NT2.Beta = Exemplar(a,b).NT2.Beta;
        else
          NE.NT2.Beta = [];
        end
        NE.NT2.Center = Exemplar(a,b).NT2.Center;
        NE.NT2.Code = Exemplar(a,b).NT2.Code;
        NE.NT2.Rot = Exemplar(a,b).NT2.Rot;
        NE.NT2.Fit = Exemplar(a,b).NT2.Fit;
        NE.NT2.Syn = Exemplar(a,b).NT2.Syn;

        NE.Count = Exemplar(a,b).Count;
        NE.Resolution = Exemplar(a,b).Resolution;
        NE.R = Exemplar(a,b).R;
        NE.T1 = Exemplar(a,b).T1;
        NE.T2 = Exemplar(a,b).T2;
        NE.AngleWeight = Exemplar(a,b).AngleWeight;
        NE.LDiscCutoff = Exemplar(a,b).LDiscCutoff;

          NE
          NE.NT1
          NE.NT2

        elseif ~isfield(Exemplar(a,b).NT1,'Unit'),

          File = zAddNTData(Exemplar(a,b).Filename,0,[],1);

          if ~isfield(File.NT(1),'Unit'),
            File = zAddNTData(Exemplar(a,b).Filename,4,[],1);
          end

          i1 = zIndexLookup(File,Exemplar(a,b).NT1.Number,Exemplar(a,b).NT1.Chain);
          i2 = zIndexLookup(File,Exemplar(a,b).NT1.Number,Exemplar(a,b).NT1.Chain);
          Exemplar(a,b).NT1 = File.NT(i1);
          Exemplar(a,b).NT2 = File.NT(i2);

File.Info.Resolution

          if isempty(File.Info.Resolution),
            Exemplar(a,b).Resolution = NaN;
          else
            Exemplar(a,b).Resolution = File.Info.Resolution;
          end

Exemplar(a,b)
Exemplar(a,b).NT1
Exemplar(a,b).NT2

        end
      end
    end
  end
