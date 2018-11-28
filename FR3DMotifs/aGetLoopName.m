%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a File structure, divides it into fragments and makes a name like
% 0/A/826:829,0/A/857:861,0/A/868:874 for a 3WJ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [loop_id] = aGetLoopName(File)

    sep = '/';
    rangesep = ':';
    fragmentsep = ',';

    br = length(File.chain_breaks);
    switch br
        case 0 %'HL'
            chbr = [1 File.NumNT];
        case 1 %'IL'
            chbr = [1 File.chain_breaks File.chain_breaks+1 File.NumNT];
        otherwise % J3, J4, etc.
            chbr = sort([1 File.chain_breaks File.chain_breaks+1 File.NumNT]);
    end

    N = length(chbr);
    fragments = cell(1,N);

    for i = 1:2:N

        fragment_start = chbr(i);
        fragment_end   = chbr(i+1);

        % NB! change this regexp when searching in multiple pdb files is implemented
%         pdbfile = regexp(File.Filename,'[a-zA-Z0-9]{4}','match');

        if isfield(File.NT(fragment_start),'ModelNum')
            model = int2str(File.NT(fragment_start).ModelNum);
        else
            model = '0'; % to distinguish from models in PDB (always start with 1)
        end

        chain = File.NT(fragment_start).Chain;
        nt1   = File.NT(fragment_start).Number;
        nt2   = File.NT(fragment_end).Number;

%         fragments{i} = [pdbfile{1} sep model sep chain sep nt1 rangesep nt2 fragmentsep];
        fragments{i} = [model sep chain sep nt1 rangesep nt2 fragmentsep];

    end

    loop_id = [fragments{:}];
    loop_id = loop_id(1:end-1); % remove the last fragmentsep

end