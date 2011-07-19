% makeAllSequencesFasta(filename,LL,RL) generates all internal loops with LL nucleotides on the left strand and RL on the right, using closing GC basepairs in every case.  

function [] = makeAllSequencesFasta(filename,LL,RL)
    L = RL + LL;
    letter = ['A','C','G','U'];
    variantNum = 4^L;
    fid = fopen(filename,'w');
    for i = 0:variantNum-1,
        fprintf(fid,'> Variant %d\n',i+1);
        fwrite(fid,'G');
        for j = 1:LL
            li = floor(i/(4^(L-j)));
            li = mod(li,4)+1;
            fprintf(fid,'%s',letter(li));
        end
        fwrite(fid,'G*C');
        for j = 1:RL
            li = floor(i/(4^(L-LL-j)));
            li = mod(li,4)+1;
            fprintf(fid,'%s',letter(li));
        end
        fprintf(fid,'C\n');
    end
    fclose(fid);
end
