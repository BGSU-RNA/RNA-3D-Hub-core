function [name, M] = aGetHetId(File, ind)

    Het = File.Het(ind);
    
    % biological assembly vs asymmetric unit
    if strfind(File.PDBFilename,'_')
        if strfind(File.PDBFilename,'_A')
            pdb_type = 'AU';
        else
            pdb_type = ['BA' File.PDBFilename(6)]; %1J5E_1.pdb
        end
    else
        if regexp(File.PDBFilename,'\d$') % .pdb1 etc
            pdb_type = ['BA' File.PDBFilename(end)];
        else
            pdb_type = 'AU';
        end
    end

    % model number. Default 1.
    if ~isfield(Het,'ModelNum')
        Het.ModelNum = 1;
    end
    M = Het.ModelNum;
    
    %pdb id
    pdb_id = regexp(File.Filename,'[a-zA-Z0-9]{4}','match');   

    % alternate location - turned off for now. Default '' or 'A'
%     if AA.AltLoc == ' ';
%         alternateId = '';
%     else
%         alternateId = NT.AltLoc;
%     end
    
    % insertion code
    insertion = regexp(Het.Number,'[a-zA-Z]','match');
    if isempty(insertion)
        insCode = '';
    else
        insCode = insertion{1};
        Het.Number = Het.Number(1:end-1);
    end

    name = sprintf('%s_%s_%i_%s_%s_%s_%s',...
    pdb_id{1},pdb_type,Het.ModelNum,Het.Chain,Het.Number,Het.Unit,insCode);
%         pdb_id{1},pdb_type,AA.ModelNum,AA.Chain,AA.Unit,AA.Number,insCode);
%     name = sprintf('%s_%s_%i_%s_%s_%s_%s_%s',pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId);
    
end