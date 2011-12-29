function [Het] = aParseHetEntities(F)

    for i = 1:length(F.Het)
        % NMR structures may not have this field
        if ~isfield(F.Het(i),'ModelNum') || isempty(F.Het(i).ModelNum)
            
            F.Het(i).ModelNum = 1;
        end
    end

    h      = F.Het(1);
    Het(1) = new_het_struct(h);
    Het    = push_info(Het, 1, h);    

    for i = 2:length(F.Het)

        h   = F.Het(i);
        ind = find_existing_het(Het, h);
        if ~isequal(ind,0)
            Het = push_info(Het, ind, h);
        else
            Het(end+1) = new_het_struct(h); %#ok<AGROW>
            Het        = push_info(Het, length(Het), h);    
        end
    
    end
    
%     keyboard;
    
end

function [H] = new_het_struct(h)

    H.Unit = h.Unit;
    H.Chain = h.Chain;
    H.Number = h.Number;
    H.Center = h.Center;
    H.ModelNum = h.ModelNum;
    
    H.Atom = {};
    H.AtomNumber = [];
    H.Beta = [];
    H.Loc = [];

end

function [Het] = push_info(Het, ind, h)

    Het(ind).Atom{end+1,1} = h.Atom;
    Het(ind).AtomNumber(end+1,1) = h.AtomNumber;
    Het(ind).Beta(end+1) = h.Beta;
    Het(ind).Loc(end+1,1:3) = h.Loc;    
    
end

function [ind] = find_existing_het(Het, h)

    ind = 0;
    
%     for i = length(Het):-1:1
        if strcmp(Het(end).Number,h.Number) && ...
           strcmp(Het(end).Chain,h.Chain)   && ...
           Het(end).ModelNum == h.ModelNum
                      
            ind = length(Het);
            return;
        end                   
%     end

end


    
    
    
%     F = zAddNTData(pdb);

%     Het = struct;
%     [a,b,c] = unique({F.Het.Number});
%     [x,y,z] = unique([F.Het.ModelNum]);
    
    
    
%     for i = 1:length(a)
%         ind = find(c == i);
%         Het(i).Unit     = F.Het(ind(1)).Unit;
%         Het(i).Chain    = F.Het(ind(1)).Chain;
%         Het(i).Number   = F.Het(ind(1)).Number;
%         Het(i).Center   = F.Het(ind(1)).Center;
%         Het(i).ModelNum = 1;
%         
%         Het(i).Atom       = cell(length(ind),1);
%         Het(i).AtomNumber = zeros(length(ind),1);
%         Het(i).Beta       = Het(i).AtomNumber;
%         Het(i).Loc        = zeros(length(ind),3);
%         
%         for j = 1:length(ind)
%             Het(i).Atom{j,1}       = F.Het(ind(j)).Atom;
%             Het(i).AtomNumber(j,1) = F.Het(ind(j)).AtomNumber;
%             Het(i).Beta(j,1)       = F.Het(ind(j)).Beta;
%             Het(i).Loc(j,1:3)      = F.Het(ind(j)).Loc;
%         end    
%     
%     end
    
    
% %     keyboard;
%     
%     
%     fid = fopen('test_mod_nts.pdb','w');
%     a = 1;
%     for i = 1:length(Het)
%         a = aWriteModNTPDB(fid,Het(i),a);        
%     end
%     fclose(fid);
%     