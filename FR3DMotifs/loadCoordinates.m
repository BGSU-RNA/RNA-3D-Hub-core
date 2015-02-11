%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input: pdb id
% output: filename with output, exit status and error message
% exit codes:
%   0 = success
%   1 = failure
%   2 = no nucleotides in pdb file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FILENAME, status, err_msg] = loadCoordinates(pdb_id)

    try
        FILENAME = fullfile(pwd, 'Coordinates.csv');
        status   = 2;
        err_msg  = '';

        F = zAddNTData(pdb_id);

        if isempty(F.NT)
            return;
        end

        fid = fopen(FILENAME,'w');

        if isfield(F,'Het') && ~isempty(F.Het)
            F.Het = aParseHetEntities(F);
        else
            F.Het = [];
        end

        N = length(F.NT);
        A = length(F.AA);
        H = length(F.Het);

        C = struct('id',  cell(1,N+A+H), ...
                   'pdb',     '', ...
                   'pdb_type','', ...
                   'model',   '', ...
                   'chain',   '', ...
                   'number',  '', ...
                   'unit',    '', ...
                   'ins_code','', ...
                   'index',   '', ...
                   'coordinates','');

        for i = 1:N
            id   = F.NT(i).ID;
            C(i) = aParseId(id,i);
            C(i).coordinates = aGetNucleotideAsPdbText(F.NT(i));
            print_line(fid, C(i));
        end

        for i = 1:A
            id         = F.NT(i).ID;
            current    = N+i;
            C(current) = aParseId(id,i);
            C(current).coordinates = aGetAsPdbText(F.AA(i),'aa');
            print_line(fid, C(current));
        end

        S = N + A;
        for i = 1:H
            id         = F.NT(i).ID;
            current    = S+i;
            C(current) = aParseId(id,i);
            C(current).coordinates = aGetAsPdbText(F.Het(i),'het');
            print_line(fid, C(current));
        end

        fclose(fid);
        status = 0;

    catch err
        err_msg = sprintf('Error "%s" in loadCoordinates on line %i (%s)\n', err.message, err.stack.line, pdb_id);
        disp(err_msg);
        status = 1;
    end

end

function [line] = print_line(fid,C)

    line = fprintf(fid,'"%s","%s","%s","%s","%s","%s","%s","%s","%i","%s"\n', ...
                   C.id, C.pdb, C.pdb_type, C.model, ...
                   C.chain, C.number, C.unit, C.ins_code, C.index, ...
                   C.coordinates);

end

function [Coord] = aParseId(id, index)

    parts          = regexp(id,'_','split');
    Coord.id       = id;
    Coord.pdb      = parts{1};
    Coord.pdb_type = parts{2};
    Coord.model    = parts{3};
    Coord.chain    = parts{4};
    Coord.number   = parts{5};
    Coord.unit     = parts{6};
    Coord.ins_code = parts{7};
    Coord.index    = index;
    Coord.coordinates = '';

end

function [Text] = aGetAsPdbText(AA, aa_or_het)

    N = length(AA.Atom(:,1));

    if strcmp(aa_or_het,'het')
        prefix = 'HETATM %4d';
    elseif strcmp(aa_or_het,'aa')
        prefix = 'ATOM %6d';
    else
        error('Incorrect mode');
    end

    lines = cell(1,N);

    for i = 1:N
        lines{i} = sprintf([prefix '  %-3s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f\n'], ...
                            i, AA.Atom{i}, AA.Unit, AA.Chain, AA.Number, ...
                            AA.Loc(i,1), AA.Loc(i,2), AA.Loc(i,3), 1, AA.Beta(i));
    end

    Text = [lines{:}];
    Text = strtrim(Text);

end

% slightly modified zWriteNucleotidePDB.m
function [Text] = aGetNucleotideAsPdbText(NT)

    Text = '';
    c = 0;
    R = eye(3);
    Sh = [0 0 0];
    a = 1;

    x = mod(c,30);                       % shift along x axis
    y = mod(fix(c/30),30);               % shift along y axis
    z = mod(fix(c/900),900);             % shift along z axis

    Lim(1,:) = [10  8 11 8  Inf Inf];   % number of base atoms, excluding hydrogen
    Lim(2,:) = [15 13 16 12 Inf Inf];   % total number of atoms, including hydrogen
    Lim(3,:) = [13  9 14 10 Inf Inf];   % locations of fictitious hydrogens

    A = {' N9' ' C4' ' N3' ' N1' ' C6' ' N6' ' C8' ' C5' ' C2' ' N7' ' H2' ' H8' ' H9' '1H6' '2H6'};
    C = {' N1' ' C2' ' O2' ' N3' ' C4' ' N4' ' C6' ' C5' ' H1' ' H6' ' H5' '1H4' '2H4'};
    G = {' N9' ' C4' ' N3' ' N1' ' C6' ' O6' ' C8' ' C5' ' C2' ' N7' ' N2' ' H1' ' H8' ' H9' '1H2' '2H2'};
    U = {' N1' ' C2' ' O2' ' N3' ' C4' ' O4' ' C6' ' C5' ' H5' ' H1' ' H3' ' H6'};
    % S = {' C1*' ' C2*' ' O2*' ' C3*' ' O3*' ' C4*' ' O4*' ' C5*' ' O5*' ' P  ' ' O1P' ' O2P'};
    S = {' C1*' ' C2*' ' O2*' ' C3*' ' O3*' ' C4*' ' O4*' ' C5*' ' O5*' ' P  ' ' OP1' ' OP2'};

    A = strrep(A,'*','''');
    C = strrep(C,'*','''');
    G = strrep(G,'*','''');
    U = strrep(U,'*','''');
    S = strrep(S,'*','''');

    %      1     2     3     4     5     6     7     8     9     10  11    12
    SugarReorder = [10 11 12 9 8 6 7 4 5 2 3 1];

    for k = 1:min(12,length(NT.Sugar(:,1))),         % loop through sugar atoms
        j = SugarReorder(k);
        L = (NT.Sugar(j,:) - Sh)*R;
        L = L + 30*[x y z];
        Text = [Text, sprintf('ATOM %6d', a)];
        Text = [Text, sprintf(' %3s', S{j})];
        Text = [Text, sprintf(' %3s', NT.Base)];
        Text = [Text, sprintf(' %1s', NT.Chain)];
        Text = [Text, sprintf('%4s    ', NT.Number)];
        Text = [Text, sprintf('%8.3f', L)];
        Text = [Text, sprintf('%6.2f', 1)];
        Text = [Text, sprintf('%6.2f\n', 99.99)];
        a = a + 1;
    end
    s = length(NT.Fit(:,1));

    for j = 1:min(Lim(2,NT.Code),s),             % loop through base atoms and H
        if j ~= Lim(3,NT.Code) && NT.Code ~= 6,
            Text = [Text, sprintf('ATOM %6d', a)];
            switch NT.Code,
                case 1, Text = [Text, sprintf(' %3s', A{j})];
                case 2, Text = [Text, sprintf(' %3s', C{j})];
                case 3, Text = [Text, sprintf(' %3s', G{j})];
                case 4, Text = [Text, sprintf(' %3s', U{j})];
                case 5, Text = [Text, sprintf(' %3s', NT.AtomName{j})];  % modified
                case 6, Text = [Text, sprintf(' %3s', NT.AtomName{j})];  % missing base atom(s)
            end
            Text = [Text, sprintf('  %3s', NT.Base)];
            Text = [Text, sprintf(' %1s',  NT.Chain)];
            Text = [Text, sprintf('%4s    ',   NT.Number)];
            L = (NT.Fit(j,:) - Sh)*R;
            L = L + 30*[x y z];                        % shift 20 Angstroms
            Text = [Text, sprintf('%8.3f', L)];        % write atom location
            Text = [Text, sprintf('%6.2f', 1)];
            Text = [Text, sprintf('%6.2f\n', 99.99)];
            a = a + 1;
        end
    end

    Text = strtrim(Text);

end

% COLUMNS        DATA  TYPE    FIELD        DEFINITION
% -------------------------------------------------------------------------------------
%  1 -  6        Record name   "ATOM  "
%  7 - 11        Integer       serial       Atom  serial number.
% 13 - 16        Atom          name         Atom name.
% 17             Character     altLoc       Alternate location indicator.
% 18 - 20        Residue name  resName      Residue name.
% 22             Character     chainID      Chain identifier.
% 23 - 26        Integer       resSeq       Residue sequence number.
% 27             AChar         iCode        Code for insertion of residues.
% 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
% 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
% 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
% 55 - 60        Real(6.2)     occupancy    Occupancy.
% 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
% 77 - 78        LString(2)    element      Element symbol, right-justified.
% 79 - 80        LString(2)    charge       Charge  on the atom.
