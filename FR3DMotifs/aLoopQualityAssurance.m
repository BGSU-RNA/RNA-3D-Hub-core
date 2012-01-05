%
% codes: 1 valid, 2 missing, 3 modified, 4 self-complementary, 5 abnormal
%

function [result, L, err_msg] = aLoopQualityAssurance(pdb_id)

    try
        result  = struct();
        err_msg = '';
        L = 0; 

        list = dir(fullfile('MotifAtlas','PrecomputedData',pdb_id,'*.mat'));
        L = length(list);        
        result = struct('id', cell(1,L), ...
                        'status', 0, ...
                        'self_complementary', 0, ...
                        'modres','');
        for i = 1:L
            load(fullfile('MotifAtlas','PrecomputedData',pdb_id,list(i).name));
    %         load(fullfile('MotifAtlas','PrecomputedData',loop_id(4:7),loop_id));
            result(i).id = list(i).name(1:end-4);
            [result(i).status,result(i).modres] = aCheckAllBreaks(File);
            result(i).self_compl = aCheckILSelfComplementarity(File);
            if result(i).status == 1 && result(i).self_compl == 1
                result(i).status = 0;
            end
        end

    catch err
        err_msg = sprintf('Error "%s" on line %i (%s)\n', err.message, err.stack.line, pdb_id);
        disp(err_msg);      
    end    

end

function [allowed] = aGetAllowedChainBreaks(name)

    switch name(1)
        case 'I' % IL
            allowed = 1;
        case 'H' % HL
            allowed = 0;
        case 'J'
            allowed = str2double(name(2)) - 1; % J3_
        otherwise
            error('Cannot identify loop type');
    end    
    
end

function [code, modres] = aCheckAllBreaks(File)

    code   = -1;
    modres = '';

    breaks         = aFindBreaks(File);    % based on .Covalent matrix     
    allowed_breaks = aGetAllowedChainBreaks(File.Filename);
    
    if length(breaks) == allowed_breaks        
        fprintf('All breaks are allowed\n');
        code = 1;        
        return;
    end

    if strcmp(File.Filename(1:2),'HL') && length(unique({File.NT.Chain})) > 1
        code = 5; %Compromised hairpin loop!
        return;
    end
    
%     flank_breaks   = aFindIntervals(File); % based on .Flank matrix
    flank_breaks = File.chain_breaks;
    breaks       = setdiff(breaks,flank_breaks);
    real_breaks  = 0;
    all_mod = {};

    fprintf('Suspicious loop,%s\n',aImplode({File.NT.Number}));

    pdbfile = aParsePdbFilename(File.Filename);    
    for i = 1:length(breaks)

        nt1 = File.NT(breaks(i));
        nt2 = File.NT(breaks(i)+1);               
        
        if nt1.Chain == nt2.Chain
            [a,b] = sort([str2double(nt1.Number) str2double(nt2.Number)]);
            if isequal(b, [2,1])
                fprintf('Reordered nucleotides\n');
                temp = nt1;
                nt1  = nt2;
                nt2  = temp;
            end        
        end
        
        [status, modres] = aCheckAdjacencyInPDB(pdbfile, nt1, nt2);
        if status == 0
            error('Nucleotide not found in pdb file');
        end        
               
        if ~isempty(modres)
            all_mod = [all_mod modres]; %#ok<AGROW>            
        end
        
        if status < 0
            real_breaks  = real_breaks + 1;
        end                                                
    end
        
    if real_breaks > allowed_breaks || ~isempty(all_mod)        
        if ~isempty(all_mod)
            modres = aImplode(all_mod);
            code = 3;
            fprintf('Modified nucleotides %s, %s\n',File.Filename,modres);
        else
            code = 2;
            fprintf('Missing nucleotides %s\n',File.Filename);
        end
    else
        code = 1;
        fprintf('Valid loop\n');
    end

end

function [self_compl] = aCheckILSelfComplementarity(File)

    self_compl = 0; % not self-complementary
    if ~strcmp(File.Filename(1:2),'IL')
        return;
    end
    
%     chbr = aFindIntervals(File);
    chbr = File.chain_breaks;
    
    seq1 = [File.NT(1:chbr).Base];
    seq2 = [File.NT(chbr+1:end).Base]; 
    
    % sequences must be of the same length
    l1 = length(seq1);
    l2 = length(seq2);
    if l1 ~= l2
        return;
    end
    
    % no true non-wc pairs
    nWC_pairs = 2:12;
    File.Edge = fix(abs(File.Edge));        
    if ~isempty(intersect(File.Edge(:),nWC_pairs))
        return;    
    end    
    
    % all opposing nucleotides match canonical pairs +GU/UG
    wc = {'AU','UA','CG','GC','GU','UG'};
    for i = 1:length(seq1)
        if ~ismember([seq1(i) seq2(end-i+1)], wc)
            return;
        end        
    end

    % all conditions met    
    fprintf('Self-complementary loop %s: %s,%s\n',File.Filename,seq1,seq2);
    self_compl = 1;

end

function [pattern] = make_pattern(NT)

    if isnan(str2double(NT.Number))
        pattern = [sprintf('%3s',NT.Base) ' ' NT.Chain sprintf('%5s',NT.Number)];        
    else
        pattern = [sprintf('%3s',NT.Base) ' ' NT.Chain sprintf('%4s',NT.Number) ' '];
    end

end

function [status, modres] = aCheckAdjacencyInPDB(PDB,NT1,NT2)

    status = 0;
    modres = {};
    T = aReadPDBFile(PDB);
    if isempty(T)
        return;
    end
    
    before = make_pattern(NT1);
    after  = make_pattern(NT2);
    
    i = 1;
    while i < length(T)
        if strcmp(T(i,1:4),'ATOM') && strcmp(T(i,18:27),before)
            while strcmp(T(i,18:27),before)
                i = i + 1;
            end
            next_line = T(i,:);
            % I = inosine, DU = deoxiribose U, DNA form            
            if strcmp(next_line(1:6),'HETATM') || ...
               strcmp(next_line(20),'I')       || ...
               strcmp(next_line(19:20),'DA')   || ...
               strcmp(next_line(19:20),'DC')   || ...
               strcmp(next_line(19:20),'DG')   || ...
               strcmp(next_line(19:20),'DU')
                
                % in order to return all modified nucleotides
                while ~strcmp(next_line(18:27), after)
                    modres{end+1} = next_line(18:27); %#ok<AGROW> % 7MG B 745
                    i = i+1;
                    next_line = T(i,:);                    
                end
                modres = unique(modres);
                modres = cellfun(@(x) x(1:3),modres,'UniformOutput',false);
                modres = strtrim(modres);
                    
                status = -1; % modified nt

            elseif strcmp(next_line(18:27),after)                
                status = 1; % next nucleotide matches NT2
            else
                status = -2; % allowed chbr or missing nt or incomplete nt
            end
            return;
        end
        i = i + 1;
    end

end

function [T] = aReadPDBFile(PDB)

    global PDBFilesLocation;
    file = [PDBFilesLocation filesep PDB '.pdb'];
    
    T = [];
    if exist(file,'file')
        try
            T = textread(file,'%80c');
        catch
            T = [];
            fid = fopen(file,'r'); % make all lines 80 characters
            if fid > 0
                L = 1;
                while L > -1
                    L = fgets(fid);
                    if L > -1
                        a = length(L);
                        if a < 80,
                            L = [L ' '*ones(1,80-a)];
                        elseif a > 80,
                            L = L(1:80);
                        end
                        T = [T; L];
                    end
                end
            end
            fclose(fid);            
        end
    else
        fprintf('Attempting to download %s from PDB\n', PDB);
        try
            c = urlread(['http://www.rcsb.org/pdb/files/' PDB '.pdb']);
            T = reshape(c',81,[])';
            T = T(:,1:80);
        catch
            fprintf('Unable to download %s from PDB\n', file);
            return;
        end
    end

end
     
% http://www.wwpdb.org/documentation/format32/sect9.html

% 1 -  6        Record name   "ATOM  "
% 7 - 11        Integer       serial       Atom  serial number.
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

