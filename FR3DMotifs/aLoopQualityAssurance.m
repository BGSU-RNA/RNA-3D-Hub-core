%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check all loops in the pdb file. Return an array of structures, its
% length and an error message.

% result.id         - loop id
% result.status     - see below
% result.modres     - empty string or a list of modified nucleotides
% result.nt_sig     - a string like: 2093,2094,2095,2096,2097,2103,2161,2162
%                     for missing or incomplete nucleotides
% result.compl      - empty or the self-complementary sequences

% status codes:
% 1 - valid loop
% 2 - missing nucleotides
% 3 - modified nucleotides
% 4 - abnormal chain number
% 5 - incomplete nucleotides
% 6 - self-complementary internal loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result, L, err_msg] = aLoopQualityAssurance(pdb_id)

    try
        % initialize return values
        result  = struct();
        err_msg = '';
        L = 0;         
        
        LOOPMATFILES = fullfile('MotifAtlas','PrecomputedData');
        
        list = dir(fullfile(LOOPMATFILES, pdb_id, '*.mat'));
        L = length(list);
        
        % preallocate memory for results
        result = struct('id', cell(1,L), 'status', 0, 'modres', '', 'nt_sig', '', 'compl', '');
                    
        % read in the pdb file as a character array
        T = readPDBFile(pdb_id);
                    
        % get an array of shortened ids of missing nucleotides
        missing_ids = readRemark465();
        
        for f = 1:L % loop over .mat files
            
            load(fullfile(LOOPMATFILES, pdb_id, list(f).name));
            
            result(f).id = File.Filename;

            covalent_breaks = aFindBreaks(File); % based on .Covalent matrix     
            flank_breaks    = File.chain_breaks; % based on the FR3D search
            % analyze only suspicious breaks
            breaks          = setdiff(covalent_breaks,flank_breaks); 
            
            % some covalent_breaks are spurious. Nts are adjacent, 
            % but the covalent matrix says that they are not.
            
            if isempty(breaks) % no suspicious breaks
                [isSelfComplementary, compl] = checkSelfComplementarity();
                if isSelfComplementary
                    result(f).status = 6;
                    result(f).compl  = compl;
                    fprintf('Self-complementary %s: %s\n',File.Filename,compl);
                else
                    result(f).status = 1;
                end
            else
                nt_sig = aImplode({File.NT.Number}); % nt signature
                fprintf('Suspicious loop %s, %s\n',File.Filename,nt_sig);                
                result(f).nt_sig = nt_sig;                    

                if foundChainAbnormalities % abnormalities in chain number
                    result(f).status = 4;
                    fprintf('Abnormal chain number\n');
                else
                    [real_breaks, modifications] = checkSuspiciousBreaks;            

                    if ~isempty(modifications)
                        modres = aImplode(modifications);
                        result(f).status  = 3;
                        result(f).modres = modres;
                        fprintf('Modified nucleotides %s\n',modres);
                    else
                        if ~isempty(missing_ids) && missingNtsFound
                            result(f).status = 2;
                            fprintf('Missing nucleotides\n');                
                        else
                            if real_breaks > 0
                                result(f).status = 5;
                                fprintf('Incomplete nucleotides\n');
                            else
                                result(f).status = 1;
                                fprintf('Valid loop\n');
                            end
                        end        
                    end
                end
            end
                                    
        end % loop over .mat files

    catch err
        err_msg = sprintf('Error "%s" on line %i (%s)\n', err.message, err.stack.line, pdb_id);
        disp(err_msg);      
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [disqualify] = missingNtsFound()

        disqualify = 0;
        for b = 1:length(breaks)

            NT = getNtsAtBreakPoint(breaks(b));            

            if checkRemark465
                disqualify = 1;
                return;
            end

        end
                
       % compare next_ids with missing ids.
        function [isMissing] = checkRemark465()

            % try to generate shortened ids of nucleotides that could follow NT
            next_ids = getNextIds();

            if isempty(intersect(missing_ids, next_ids))
                isMissing = 0;
            else
                isMissing = 1;        
            end

            % return a list of plausible next ids given any nt
            function [next_ids] = getNextIds()

                next_ids = {'', '', ''};        
                number = str2double(NT.Number);

                if isnan(number) % insertion code present
                    % increment insertion code
                    % model, chain, number, ins_code++
                    ins_code = NT.Number(end-1:end) - 0; % last character -> ascii
                    incremented_ins_code = char(ins_code+1);
                    next_ids{1} = sprintf('%i_%s_%i%s',NT.ModelNum, NT.Chain, number, incremented_ins_code);
                else
                    % increment number
                    % model, chain, number+1, ins_code = ''
                    next_ids{2} = sprintf('%i_%s_%i%s',NT.ModelNum, NT.Chain, number+1, '');

                    % same number, but with an insertion code
                    % model, chain, number, ins_code = 'a'
                    next_ids{3} = sprintf('%i_%s_%i%s',NT.ModelNum, NT.Chain, number, 'a');
                end

            end            

        end
                        
    end % missingNtsPresent

    % Return the number of real breaks + all modifications, if any.
    function [real_breaks, modifications] = checkSuspiciousBreaks()

        real_breaks = 0;
        modifications = {};

        for b = 1:length(breaks)

            [nt1, nt2] = getNtsAtBreakPoint(breaks(b));            

            [isAdjacent, modif] = checkAdjacencyInPDB(nt1, nt2);

            if ~isempty(modif)
                modifications = [modifications modif]; %#ok<AGROW>            
            end

            if ~isAdjacent
                real_breaks  = real_breaks + 1;
            end        
        end        
        
        % Check if two nucleotides, NT1 and NT2, are adjacent in the PDB text file.
        % Return codes: 
        % 1 for adjecent (need to check that there are no missing nts in between)
        % 0 for non-adjacent nts, which could be due to an allowed chain break
        % or a missing nt or an incomplete nt.
        function [isAdjacent, modres] = checkAdjacencyInPDB(NT1,NT2)

            isAdjacent = -1;
            modres = {};
        %     T = aReadPDBFile(PDB);
            if isempty(T)
                error('PDB file could not be read');
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
                            i = i + 1;
                            next_line = T(i,:);                    
                        end
                        modres = unique(modres);
                        modres = cellfun(@(x) x(1:3),modres,'UniformOutput',false);
                        modres = strtrim(modres);

                        isAdjacent = 0; % modified nt

                    elseif strcmp(next_line(18:27),after)                
                        isAdjacent = 1; % next nucleotide matches NT2
                    else
                        isAdjacent = 0; % NT1 and NT2 are not adjacent. 
                    end
                    return;
                end
                i = i + 1;
            end

            assert(status ~= -1, 'Nucleotide not found in pdb file');

            function [pattern] = make_pattern(NT)
                if isnan(str2double(NT.Number))
                    pattern = [sprintf('%3s',NT.Base) ' ' NT.Chain sprintf('%5s',NT.Number)];        
                else
                    pattern = [sprintf('%3s',NT.Base) ' ' NT.Chain sprintf('%4s',NT.Number) ' '];
                end
            end

        end        
        
    end % checkSuspiciousBreaks

    % get nucleotides before and after the chainbreak
    function [nt1, nt2] = getNtsAtBreakPoint(breakpoint)

        nt1 = File.NT(breakpoint);
        nt2 = File.NT(breakpoint+1);               

        % make sure that the order is correct
        if nt1.Chain == nt2.Chain
            [a,b] = sort([str2double(nt1.Number) str2double(nt2.Number)]);
            if isequal(b, [2,1])
                fprintf('Reordered nucleotides\n');
                temp = nt1;
                nt1  = nt2;
                nt2  = temp;
            end        
        end        

    end % getNtsAtBreakPoint

    % make sure that chain numbers make sense. For example, a hairpin can
    % not contain nucleotides from more than 1 chain, an internal loop -
    % from more than 2, a 3wj - from more than 3.
    function [isAbnormal] = foundChainAbnormalities()

        isAbnormal    = 0;
        loop_type     = File.Filename(1:2);        
        unique_chains = length(unique({File.NT.Chain}));

        if strcmp(loop_type,'HL') && unique_chains > 1
            isAbnormal = 1;
        end
        if strcmp(loop_type,'IL') && unique_chains > 2
            isAbnormal = 1;
        end
        if strcmp(loop_type,'J3') && unique_chains > 3
            isAbnormal = 1;
        end

    end % foundChainAbnormalities

    % Only for ILs. Return 0 is not self-complementary, 1 otherwise.
    function [isSelfCompl, seqs] = checkSelfComplementarity()

        isSelfCompl = 0; % by default not self-complementary
        seqs = ''; % self-cpmplementary sequences

        if ~strcmp(File.Filename(1:2),'IL')
            return;
        end

        chbr = File.chain_breaks;

        seq1 = [File.NT(1:chbr).Base];
        seq2 = [File.NT(chbr+1:end).Base]; 

        % sequences must be of the same length
        if length(seq1) ~= length(seq2)
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
        isSelfCompl = 1;
        seqs = [seq1 ',' seq2];

    end

    % return a list of shortened ids of missing nucleotides from remark 465
    function [ids] = readRemark465()
                
        ids = {};        

        if isempty(T)
            error('PDB file could not be read');
        end

        i = 1;
        while i < length(T)
            if strcmp(T(i,1:10),'REMARK 465') % found the first line of 465
                j = i + 7; % skip 7 lines with general info
                while strcmp(T(j,1:10),'REMARK 465')
                    s = T(j,1:27);                    
                    model   = s(14);
                    chain   = s(20);
                  % unit    = s(17:18);
                    number  = strtrim(s(22:26));
                    inscode = s(27);
                    
                    if model == ' '
                        model = '1';
                    end                    
                    if inscode == ' '
                        inscode = '';
                    end
                    ids{end+1} = [model '_' chain '_' number inscode]; %#ok<AGROW>
                    j = j + 1;
                end                
                return; % reached the end of remark 465, return
            end
            
            if strcmp(T(i,1:4),'ATOM') % no remark 465
                return;
            end
            i = i + 1;
        end
        
    end % readRemark465

end % LoopQualityAssurance

% return a character array with one line per row
function [T] = readPDBFile(PDB)

    file = [PDB '.pdb'];

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



% _________________________________________________________________________

% http://www.wwpdb.org/documentation/format23/remarks2.html#REMARK470

% REMARK 465                                                                      
% REMARK 465 MISSING RESIDUES                                                     
% REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE                       
% REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN               
% REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)                
% REMARK 465                                                                      
% REMARK 465   M RES C SSSEQI 
% _________________________________________________________________________

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
