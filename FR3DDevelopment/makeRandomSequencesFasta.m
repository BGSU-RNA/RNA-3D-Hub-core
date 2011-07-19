% makeRandomSequencesFasta(filename,LL,RL,N) generates N internal loops with LL nucleotides on the left strand and RL on the right, using closing GC basepairs in every case.  

function [] = makeRandomSequencesFasta(filename,LL,RL,N)
    letter = ['A','C','G','U'];
    fid = fopen(filename,'w');
    for i = 1:N,
        fprintf(fid,'> Variant %d\n',i);
        fwrite(fid,'G');
        m = ceil(4*rand(1,LL));               % uniform distribution ACGU
        fprintf(fid,'%s',letter(m));
        fwrite(fid,'G*C');
        m = ceil(4*rand(1,RL));               % uniform distribution ACGU
        fprintf(fid,'%s',letter(m));
        fprintf(fid,'C\n');
    end
    fclose(fid);
end
