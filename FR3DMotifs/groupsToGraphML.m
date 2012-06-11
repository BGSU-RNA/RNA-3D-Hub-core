function [] = groupsToGraphML(Location, groups, M, names, cutoff)

    N = length(groups);
    fid = fopen(fullfile(Location, 'Supergroups.graphml'),'w');
   
    fprintf(fid,'<graphml>\n');
    fprintf(fid,'<key id="label" for="node" attr.name="label" attr.type="string"/>\n'); %node size   
    fprintf(fid,'<key id="within" for="node" attr.name="within" attr.type="double"/>\n'); %node color    
    fprintf(fid,'<key id="cands" for="node" attr.name="cands" attr.type="double"/>\n'); %node size
    fprintf(fid,'<key id="numnt" for="node" attr.name="numnt" attr.type="double"/>\n'); %number of nucleotides
    fprintf(fid,'<key id="signature" for="node" attr.name="signature" attr.type="string"/>\n'); %motif signature
    fprintf(fid,'<key id="links" for="edge" attr.name="links" attr.type="double"/>\n'); %width
    fprintf(fid,'<key id="disc" for="edge" attr.name="disc" attr.type="double"/>\n'); %edge color    
    fprintf(fid,'<key id="connection" for="edge" attr.name="connection" attr.type="string"/>\n'); %closest loops
    fprintf(fid,'<graph edgedefault="undirected">\n');
    
%     reasons = {'basepair mismatch ','no match ','extra stacking ','extra basepair '};
    
    for i = 1:N

        load( fullfile(Location, 'Groups', sprintf('Group_%03d',i)) );        
        
        [a,b,c] = intersect(groups{i},names);        
        submatrix = M(c,c);
        
        fprintf(fid,'<node id="Group_%03d">\n',i);
        fprintf(fid,'\t<data key="label">Group_%03d</data>\n',i);
        fprintf(fid,'\t<data key="signature">%s</data>\n',Search.Signature);        
        fprintf(fid,'\t<data key="numnt">%f</data>\n',length(Search.Candidates(1,1:end-1)));                
        fprintf(fid,'\t<data key="cands">%i</data>\n',length(groups{i}));        
        fprintf(fid,'\t<data key="within">%f</data>\n',mean(submatrix(:)));
        fprintf(fid,'</node>\n');
                
        for j = (i+1):N
                        
            if i ~= j
                [x,y,z] = intersect(groups{j},names);                
                if ~isempty(find((M(c,z)<=cutoff),1))
                    submatrix = M(c,z);
                    fprintf(fid,'<edge id="%ito%i" source="Group_%03d" target="Group_%03d">\n',i,j,i,j);
                    fprintf(fid,'\t<data key="links">%i</data>\n',round(length(find(M(c,z)<=cutoff))/2));
                    
                    [minValue,ind] = min(submatrix(:));
                    [loop1,loop2] = ind2sub([length(c) length(z)],ind);
                    loop1 = names{c(loop1)};
                    loop2 = names{z(loop2)};
                    fprintf(fid,'\t<data key="connection">%s %s</data>\n',loop1,loop2);
                    
%                     minMatchMatrix = M([c(loop1) z(loop2)],[c(loop1) z(loop2)]);
%                     if ~isempty(find(minMatchMatrix==5,1))
%                         
%                     end
                    
                    fprintf(fid,'\t<data key="disc">%f</data>\n',min(submatrix(:)));
                    fprintf(fid,'</edge>\n');
                end
            end
        end
    end
    
    fprintf(fid,'</graph>');
    fprintf(fid,'</graphml>');
    fclose(fid);

end