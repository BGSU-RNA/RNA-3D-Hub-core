% zEvaluateTriples loads a SearchSaveFile, keeps one instance of each triple, makes two models for each triple (using the observed interaction between 1 and 2 and 2 and 3, ignoring any interaction between 1 and 3) and calculates discrepancies, then sorts by decreasing geometric discrepancy and displays the triples


%load 2010-07-11_00_32_13-Coplanar_triples_all_cp_ncp_excluding_regular_triples_excluding_1QCU.mat                         % not sure what this is

if 0 > 1,                                 % skip the first part for now

clear

% ----------------------------------- Part 0 - Load basic data

load 2010-07-20_11_45_16-cp_ncp_not_stack_triples_excl_1QCU

Query = Search.Query;
% File  = Search.File;

[File,FIndex] = zAddNTData(Search.Filenames,0,[],1);
File = File(FIndex);

File = rmfield(File,'AA');

VP.Sugar    = 1;
VP.AtOrigin = 2;

C = Search.Candidates;                      % candidates, with file numbers

C = double(C);                            % convert to double

for n = 1:length(C(:,1)),
  Indices = C(n,1:3);                     % indices of this triple
  f       = C(n,4);                       % file number
    
  C(n,8)  = File(f).Edge(Indices(1),Indices(2));  % observed BPs
  C(n,9)  = File(f).Edge(Indices(2),Indices(3));
  C(n,10) = File(f).Edge(Indices(3),Indices(1));

  if C(n,8)  == 0, C(n,8)  = 1000; end
  if C(n,9)  == 0, C(n,9)  = 1000; end
  if C(n,10) == 0, C(n,10) = 1000; end

  C(n,5:7) = abs(C(n,8:10));              % better for sorting
end

% round(C(1:20,:))

D = C;

D(:,1:3) = sort(C(:,1:3),2);                % sort the indices of each triple

[U,n,i] = zUniqueRows(D,1:4);               % find unique rows using cols 1:4
C = C(i,:);                                 % keep one instance of each triple

% round(C(1:20,:))

s = length(C(:,1));                         % number of triples remaining

Geo = Inf * ones(s,1);                      % place to store geom discreps
Iso = Inf * ones(s,1);

VPP.AtOrigin = 2;
VPP.LabelBases = 0;
VPP.Sugar = 0;
VPP.LineStyle = '-.';
VPP.Title = 0;

% ----------------------------- Part I - Make and score models for each triple


for v = 1:2,                                % pass one, compute discreps
                                            % pass two, display triples

  for n = 1:s,

    fprintf('\nTriple %d of %d\n',n,s);

    Indices = C(n,1:3);                     % indices of this triple
    f       = C(n,4);                       % file number

    N1 = File(f).NT(Indices(1));
    N2 = File(f).NT(Indices(2));
    N3 = File(f).NT(Indices(3));

    BP12 = File(f).Edge(Indices(1),Indices(2));  % observed BP interaction
    BP23 = File(f).Edge(Indices(2),Indices(3));
    BP31 = File(f).Edge(Indices(3),Indices(1));

    Nam = ['_Res_' strrep(num2str(File(f).Info.Resolution),'.',',') '_' File(f).Filename '_' N1.Base N1.Number '_' N2.Base N2.Number '_' N3.Base N3.Number '_'];

    Nam = [Nam zEdgeText(BP12,2) '_' zEdgeText(-BP31,2) '_' zEdgeText(BP23,2) '_Model_'];

    % convert basepair classifications from near to true

    if BP12 > 100,  BP12 = BP12 - 100; end
    if BP12 < -100, BP12 = BP12 + 100; end
    if BP23 > 100,  BP23 = BP23 - 100; end
    if BP23 < -100, BP23 = BP23 + 100; end
    if BP31 > 100,  BP31 = BP31 - 100; end
    if BP31 < -100, BP31 = BP31 + 100; end

    if exist('comparetoregular.txt') > 0,     % don't use subcategories
      BP12 = fix(BP12);
      BP23 = fix(BP23);
      BP31 = fix(BP31);
    end

    Nam1213 = [Nam zEdgeText(BP12,2) '_' zEdgeText(-BP31,2) '_----_'];
    Nam1223 = [Nam zEdgeText(BP12,2) '_----_' zEdgeText(BP23,2) '_'];

    Nam1223 = strrep(Nam1223,' ','');
    Nam1213 = strrep(Nam1213,' ','');

    if v > 1,                             % display triples on the second pass
      figure(10)
      clf
      cla
      figure(11)
      clf
      cla
    end

    if BP12 ~= 0 && BP23 ~= 0 && BP31 ~= 0,
      FFF = zMakeTriple(BP12,BP23,N1.Code,N2.Code,N3.Code);
      if ~isempty(FFF),
        hold on
        if v > 1,
          figure(10)
          clf
          cla
          zDisplayNT(FFF,1:3,VPP)
          ytext = 'Model interactions ';
          ytext = [ytext ' 12:' zEdgeText(BP12,1) ' '];
          ytext = [ytext ' 23:' zEdgeText(BP23,1) ' '];
%         ytext = [ytext ' 13:' zEdgeText(BP31,1) ' '];

          zDisplayNT(File(f),Indices,VP);
          xlabel(ytext);
          view(2)
        end

        fprintf('12-23\n');
        fprintf('%s\n',Nam1223);
        GeoD = xDiscrepancy(File(f),Indices(1:3),FFF,1:3);
        IsoD = xDiscrepancyForTriples(File(f),Indices(1:3),FFF,1:3);
        Geo(n) = min(Geo(n),GeoD);
        Iso(n) = min(Iso(n),IsoD);
        fprintf('Geometric discrepancy from model: %7.4f\n',GeoD);
        fprintf('Triple IsoDiscrepancy from model: %7.4f\n',IsoD);

        if v > 1 && GeoD == Geo(n),
          G = strrep(num2str(GeoD),'.',',');
          rotate3d on
          saveas(gcf,['Triples' filesep G Nam1223 '.fig'],'fig');
          saveas(gcf,['Triples' filesep G Nam1223 '.png'],'png');
        end
      end

      FFF = zMakeTriple(BP31,BP12,N3.Code,N1.Code,N2.Code);
      if ~isempty(FFF),
        hold on
        FFF.NT = FFF.NT([2 3 1]);

        if v > 1,
          figure(11)
          clf
          cla
          zDisplayNT(FFF,1:3,VPP)
          ytext = 'Model interactions ';
          ytext = [ytext ' 12:' zEdgeText(BP12,1) ' '];
%         ytext = [ytext ' 23:' zEdgeText(BP23,1) ' '];
          ytext = [ytext ' 13:' zEdgeText(BP31,1) ' '];
          zDisplayNT(File(f),Indices,VP);
          xlabel(ytext);
          view(2)
          axis vis3d
          axis equal
        end

        fprintf('12-13\n');
        fprintf('%s\n',Nam1213);
        GeoD = xDiscrepancy(File(f),Indices(1:3),FFF,1:3);
        IsoD = xDiscrepancyForTriples(File(f),Indices(1:3),FFF,1:3);
        Geo(n) = min(Geo(n),GeoD);
        Iso(n) = min(Iso(n),IsoD);
        fprintf('Geometric discrepancy from model: %7.4f\n',GeoD);
        fprintf('Triple IsoDiscrepancy from model: %7.4f\n',IsoD);
        if v > 1 && GeoD == Geo(n),
          G = strrep(num2str(GeoD),'.',',');
          rotate3d on
          saveas(gcf,['Triples' filesep G Nam1213 '.fig'],'fig');
          saveas(gcf,['Triples' filesep G Nam1213 '.png'],'png');
        end
      end

    elseif BP12 ~= 0 && BP23 ~= 0,
      FFF = zMakeTriple(BP12,BP23,N1.Code,N2.Code,N3.Code);
      if ~isempty(FFF),
        hold on
        if v > 1,
          figure(10)
          clf
          cla
          zDisplayNT(FFF,1:3,VPP)
          ytext = 'Model interactions ';
          ytext = [ytext ' 12:' zEdgeText(BP12,1) ' '];
          ytext = [ytext ' 23:' zEdgeText(BP23,1) ' '];
%         ytext = [ytext ' 13:' zEdgeText(BP31,1) ' '];
          xlabel(ytext);

          zDisplayNT(File(f),Indices,VP);
          view(2)
        end

        fprintf('12-23\n');
        fprintf('%s\n',Nam1223);
        GeoD = xDiscrepancy(File(f),Indices(1:3),FFF,1:3);
        IsoD = xDiscrepancyForTriples(File(f),Indices(1:3),FFF,1:3);
        Geo(n) = min(Geo(n),GeoD);
        Iso(n) = min(Iso(n),IsoD);
        fprintf('Geometric discrepancy from model: %7.4f\n',GeoD);
        fprintf('Triple IsoDiscrepancy from model: %7.4f\n',IsoD);
        if v > 1 && GeoD == Geo(n),
          G = strrep(num2str(GeoD),'.',',');
          rotate3d on
          saveas(gcf,['Triples' filesep G Nam1223 '.fig'],'fig');
          saveas(gcf,['Triples' filesep G Nam1223 '.png'],'png');
        end
      end

    elseif BP12 ~= 0 && BP31 ~= 0,
      FFF = zMakeTriple(BP31,BP12,N3.Code,N1.Code,N2.Code);
      if ~isempty(FFF),
        hold on
        FFF.NT = FFF.NT([2 3 1]);
         if v > 1,
          figure(11)
          clf
          cla
          zDisplayNT(FFF,1:3,VPP)
          ytext = 'Model interactions ';
          ytext = [ytext ' 12:' zEdgeText(File(f).Edge(Indices(1),Indices(2)),1) ' '];
%         ytext = [ytext ' 23:' zEdgeText(File(f).Edge(Indices(2),Indices(3)),1) ' '];
          ytext = [ytext ' 13:' zEdgeText(File(f).Edge(Indices(1),Indices(3)),1) ' '];
          xlabel(ytext);
          zDisplayNT(File(f),Indices,VP);
          view(2)

        end

        fprintf('12-13\n');
        fprintf('%s\n',Nam1213);
        GeoD = xDiscrepancy(File(f),Indices(1:3),FFF,1:3);
        IsoD = xDiscrepancyForTriples(File(f),Indices(1:3),FFF,1:3);
        Geo(n) = min(Geo(n),GeoD);
        Iso(n) = min(Iso(n),IsoD);
        fprintf('Geometric discrepancy from model: %7.4f\n',GeoD);
        fprintf('Triple IsoDiscrepancy from model: %7.4f\n',IsoD);
        if v > 1 && GeoD == Geo(n),
          G = strrep(num2str(GeoD),'.',',');
          rotate3d on
          saveas(gcf,['Triples' filesep G Nam1213 '.fig'],'fig');
          saveas(gcf,['Triples' filesep G Nam1213 '.png'],'png');
        end
      end
    else
        hold on
        if v > 1,
          figure(10)
          clf
          cla
          zDisplayNT(File(f),Indices,VP);
          view(2)
        end

        fprintf('12-23\n');
        fprintf('%s\n',Nam1223);
        if v > 1,
          G = '999';
          rotate3d on
          saveas(gcf,['Triples' filesep G Nam1223 '.fig'],'fig');
          saveas(gcf,['Triples' filesep G Nam1223 '.png'],'png');
        end
    end

    if v > 1,
      F.NT = File(f).NT(Indices);
      NT = F.NT(2);
      zWritePDB(F,['Triples' filesep G Nam1223 '.pdb'],NT.Rot,NT.Fit(1,:));
    end

    if v > 1,
%      pause
    end

    if Geo(n) == Inf,
      drawnow
      Nam
      [BP12 -BP31 BP23]
%      pause
    end
  end

  CC = C;

  if v == 1,
    [y,i] = sort(-Geo);
    C   = C(i,:);
    Geo = Geo(i);
    Iso = Iso(i);

CCC = C;

    C = [C Geo Iso];

CCCC = C;

    save TripleCoplanarData C

CCCCC = C;

  end

  fprintf('===========================================\n');

  % Display the whole list, sorted by decreasing geom. discrepancy

  for n = 1:s,
    Indices = C(n,1:3);                     % indices of this triple
    f       = C(n,4);                       % file number

    N1 = File(f).NT(Indices(1));
    N2 = File(f).NT(Indices(2));
    N3 = File(f).NT(Indices(3));

    BP12 = File(f).Edge(Indices(1),Indices(2));  % observed BP interaction
    BP23 = File(f).Edge(Indices(2),Indices(3));
    BP13 = File(f).Edge(Indices(1),Indices(3));
    
    fprintf('Triple %4d geom. discrepancy %7.4f IDI %7.4f %4s %s%4s-%s%4s-%s%4s %4s %4s %4s\n', n, Geo(n), Iso(n), File(f).Filename, N1.Base, N1.Number, N2.Base, N2.Number, N3.Base, N3.Number, zEdgeText(BP12,2), zEdgeText(BP13,2), zEdgeText(BP23,2));
  end

end

end

% ----------------------------------- Part II - FR3D searches for each triple

if 0 > 1,
  load TripleCoplanarData
  CC = C;
elseif 10 > 2,
  load TripleCoplanarDataSearch
end

  s = length(CC(:,1));                         % number of triples remaining

if exist('File'),
  File = zAddNTData('Nonredundant_4A_2010-05-19_list.pdb',0,File,1);
else
  File = zAddNTData('Nonredundant_4A_2010-05-19_list.pdb',0,[],1);
end

  load 2010-07-20_11_45_16-cp_ncp_not_stack_triples_excl_1QCU

  SS = Search;

  Verbose = 0;

  [File,EIndex] = zAddNTData(SS.Filenames,0,File,Verbose); %add PDB data if needed

diary TripleGeometricSearches.txt

%  for n = 1:s,
%  for n = 44:s,                              % edit this line for where to start

while 10 > 1,

n = input('Which triple should we look at? ');

    Indices = CC(n,1:3);                     % indices of this triple
    f       = CC(n,4);                       % file number

    N1 = SS.File(f).NT(Indices(1));
    N2 = SS.File(f).NT(Indices(2));
    N3 = SS.File(f).NT(Indices(3));

    BP12 = SS.File(f).Edge(Indices(1),Indices(2));  % observed BP interaction
    BP23 = SS.File(f).Edge(Indices(2),Indices(3));
    BP31 = SS.File(f).Edge(Indices(3),Indices(1));
    BP13 = SS.File(f).Edge(Indices(1),Indices(3));

    a = (BP12 > -13 && BP12 < 13 && BP12 ~= 0) + (BP23 > -13 && BP23 < 13 && BP23 ~= 0) + (BP13 > -13 && BP13 < 13 && BP13 ~= 0);

    if a < 200,                                % only those with < 2 basepairs

    Nam = ['_Res_' strrep(num2str(SS.File(f).Info.Resolution),'.',',') '_' SS.File(f).Filename '_' N1.Base N1.Number '_' N2.Base N2.Number '_' N3.Base N3.Number '_'];

    Nam = [Nam zEdgeText(BP12,2) '_' zEdgeText(-BP31,2) '_' zEdgeText(BP23,2) '_Model_'];

    % --------- set up a FR3D query for this triple

    clear Query
    Query.Description  = ['Triple ' Nam];
    Query.Filename     = SS.File(f).Filename;
    Query.NTList       = {N1.Number, N2.Number, N3.Number};
    Query.ChainList    = {N1.Chain,  N2.Chain,  N3.Chain};
    Query.Mask         = [N1.Base N2.Base N3.Base];
    Query.Edges{1,2}   = '~stack';
    Query.Edges{2,3}   = '~stack';
    Query.Edges{1,3}   = '~stack';
    Query.DiscCutoff   = 0.5;
    Query.SearchFiles  = 'Nonredundant_4A_2010-05-19_list';
    Query.ExcludeOverlap = 1;
    Query = xConstructQuery(Query,File(EIndex(f)));
    UsingLibrary = 1;
    Filenames = 'Nonredundant_4A_2010-05-19_list';

    xFR3DSearch

    CC(n,13) = sum(Search.Discrepancy <= 0.05);
    CC(n,14) = sum(Search.Discrepancy <= 0.1);
    CC(n,15) = sum(Search.Discrepancy <= 0.15);
    CC(n,16) = sum(Search.Discrepancy <= 0.2);
    CC(n,17) = sum(Search.Discrepancy <= 0.25);
    CC(n,18) = sum(Search.Discrepancy <= 0.3);
    CC(n,19) = sum(Search.Discrepancy <= 0.35);
    CC(n,20) = sum(Search.Discrepancy <= 0.4);
    CC(n,21) = sum(Search.Discrepancy <= 0.45);
    CC(n,22) = sum(Search.Discrepancy <= 0.5);

%n
%Nam
%round(CC(n,13:22))

    if mod(n,100) == 0,
      save TripleCoplanarDataSearchUpdate CC
    end

fprintf('Triple %d ==========================================================\n', n);

    f       = CC(n,4);                       % file number

      fprintf('Triple %4d geom. discrepancy %7.4f IDI %7.4f %4s %s%5s-%s%5s-%s%5s %5s %5s %5s %3d %3d %3d\n', n, CC(n,11), CC(n,12), SS.File(f).Filename, N1.Base, N1.Number, N2.Base, N2.Number, N3.Base, N3.Number, zEdgeText(BP12,2), zEdgeText(BP13,2), zEdgeText(BP23,2));

    xListCandidates(Search);

    end

    if 10 > 1,
      figure(1)
      clf
      cla
      rotate3d on
      xDisplayCandidates(File(SIndex),Search);
    end
  end

diary off


save TripleCoplanarDataSearchUpdate CC

% ----------------------------------- Part III - list the triples

load TripleCoplanarDataSearch
load 2010-07-20_11_45_16-cp_ncp_not_stack_triples_excl_1QCU

File = zAddNTData('Nonredundant_4A_2010-05-19_list.pdb',0,[],1)
[File,FIndex] = zAddNTData(SS.Filenames,0,File,1);  % why does it have more?

Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

  SS = Search;

s = length(CC(:,1));                         % number of triples remaining


fprintf('Triples that are not made of two basepairs\n\n');

j = [];

  for n = 1:s,
    Indices = CC(n,1:3);                     % indices of this triple
    f       = CC(n,4);                       % file number

    N1 = SS.File(f).NT(Indices(1));
    N2 = SS.File(f).NT(Indices(2));
    N3 = SS.File(f).NT(Indices(3));

    BP12 = SS.File(f).Edge(Indices(1),Indices(2));  % observed BP interaction
    BP23 = SS.File(f).Edge(Indices(2),Indices(3));
    BP13 = SS.File(f).Edge(Indices(1),Indices(3));

    for d = 1:3,
      C    = File(f).NT(Indices(d)).Code;
      b(d) = round(median(File(f).NT(Indices(d)).Beta(1:Lim(1,C))));
      CC(n,22+d) = b(d);
    end

    a = (BP12 > -13 && BP12 < 13 && BP12 ~= 0) + (BP23 > -13 && BP23 < 13 && BP23 ~= 0) + (BP13 > -13 && BP13 < 13 && BP13 ~= 0);

    if a < 2,

      if ~strcmp(SS.File(f).Filename,File(f).Filename),
        fprintf('*************** File name mismatch %d %s %s\n', f, SS.File(f).Filename, File(f).Filename);
      end

      j = [j n];                      % accumulate such triples

      fprintf('Triple %4d geom. discrepancy %7.4f IDI %7.4f %4s %s%5s-%s%5s-%s%5s %5s %5s %5s %3d %3d %3d', n, CC(n,11), CC(n,12), SS.File(f).Filename, N1.Base, N1.Number, N2.Base, N2.Number, N3.Base, N3.Number, zEdgeText(BP12,2), zEdgeText(BP13,2), zEdgeText(BP23,2), b(1), b(2), b(3));

      for k = 13:22,
        fprintf('%4d ', CC(n,k));
      end
      fprintf('\n');
    end
  end

figure(5)
clf
%plot(CC(:,23),CC(:,19),'.');
plot(CC(j,23),CC(j,19),'.');
xlabel('Median b-factor of NT1');
ylabel('Number of instances of this triple found');
title('Triples with fewer than two basepairs');
saveas(gcf,'Instances versus b-factor','png')
