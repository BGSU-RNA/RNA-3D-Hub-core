%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over pdb files downloaded from RLooM webserver and match them to 
% their locations in 1S72 using geometric FR3D searches. This procedure is
% necessary because RLooM uses an internal indexing system when referring 
% to nucleotides in the loops.

% aRloomLoops('path_to/loop_extraction_benchmark/rloom/structures')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = aRloomLoops(Location)

    pdblist = '1S72';  
    Disc = 0.001;
    fid = fopen([Location filesep 'loops.txt'],'w');

    list = dir([Location filesep '*.pdb']);
    
    F = zAddNTData(pdblist);
    
    for f = 1:length(list)
        clear Search Query;
        matfile = [Location filesep list(f).name(1:end-4) '.mat'];
        searchfile = [Location filesep list(f).name(1:end-4) '_search.mat'];
        
        if ~exist(searchfile, 'file')
            [Query, G] = constructGeometricQuery([Location filesep list(f).name]);            
%             if Query.NumNT < 25
disp(Query.NumNT);
disp(Query.NTList);
if Query.NumNT > 50
    fprintf('%s t oo big, %i nucleotides\n', Query.Name, Query.NumNT);
    continue;
end
                Search = aFR3DSearch(Query, G);                
                if isfield(Search,'Candidates') && ~isempty(Search.Candidates)
                    try
                    loop_id = makeLoopName();
                    catch,keyboard;end
                    disp(loop_id);
                    fprintf(fid, '"%s","%s"\n', loop_id, chain);        
                else
                    save(searchfile, 'Search');
                    fprintf('%s no candidates\n', Search.Query.Name);            
                end                
%             else
%                 fprintf(fid, '%s too big\n', Query.Name);
%             end
        else
            load(searchfile);
            [loop_id, chain] = makeLoopName();
            if ~strcmp(chain, '')
                fprintf(fid, '"%s","%s"\n', loop_id, chain);        
            end
            disp(loop_id)
        end
    end
    
    fclose(fid);

    function [loop_id, chain] = makeLoopName()


        
        if ~isfield(Search,'Candidates')
            loop_id = 'No results';
            chain = '';
            return;
        end
        chain = Search.File.NT(Search.Candidates(1,1)).Chain;
        loop_id = aImplode({Search.File.NT(Search.Candidates(1,1:end-1)).Number});
%         
%         sep = '/';
%         rangesep = ':';
%         fragmentsep = ',';
% 
%         pdbfile = Search.Candidates(1,end);
%         nts = Search.Candidates(1,1:end-1);
%         
%         chain_breaks = find(diff(Search.Candidates(1:end-1))>1);        
%         
%         File = Search.File;
%         br = length(chain_breaks);
%         switch br
%             case 1 %'IL'
%                 chbr = [nts(1) nts(chain_breaks) nts(chain_breaks+1) nts(end)];
%             case 0 %'HL'
%                 chbr = [nts(1) nts(end)];
%             case 2 %'J3'
%                 chbr = [nts(1) nts(chain_breaks(1)) ...
%                         nts(chain_breaks(1)+1) nts(chain_breaks(2)) ...
%                         nts(chain_breaks(2)+1) nts(end)];                
%             otherwise
%                 loop_id = 'Not a recognized loop motif';
%                 return
%         end
% 
%         N = length(chbr);
%         fragments = cell(1,N);
%         model = '1';
%         
%         for i = 1:2:N
% 
%             fragment_start = chbr(i);        
%             fragment_end   = chbr(i+1);
% 
%             chain = File.NT(fragment_start).Chain;
%             nt1   = File.NT(fragment_start).Number;
%             nt2   = File.NT(fragment_end).Number;        
% 
%             fragments{i} = [model sep chain sep nt1 rangesep nt2 fragmentsep];                
% 
%         end
% 
%         loop_id = [fragments{:}];
%         loop_id = loop_id(1:end-1); % remove the last fragmentsep                
    end    
    
    function [Query, G] = constructGeometricQuery(file)

        Query.DiscCutoff  = Disc;
        Query.SearchFiles = pdblist;    
        Query.SaveDir     = Location;      
        Query.Filenames   = pdblist;

        if exist(matfile,'file')
            load(matfile);
        else
            [File,QIndex] = zAddNTData(file(1:end-4));
            save(matfile, 'File');
        end
        G.File = [File F]; % 1 - query, 2 - target
        G.QIndex = [1 2];        
                
        m = regexp(File.Filename, filesep, 'split');
        File.Filename = m{end};
        
        Query.Name      = [File.Filename '_search.mat'];    
        Query.Geometric = 1;
        Query.NumNT     = length(File.NT);
        Query.ExcludeOverlap = 1;

        Query.ChainList = {File.NT.Chain};
        Query.NTList = {File.NT.Number};
        
        Query.Diff  = cell(Query.NumNT,Query.NumNT - 1);
        Query.Edges = cell(Query.NumNT,Query.NumNT);

    end

end