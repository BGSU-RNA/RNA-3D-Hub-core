function [] = makeSequencesFasta(filename,LL,RL,num)
    L = RL + LL;
    letter = ['A','C','G','U'];
    varientNum = 4^L;
    if nargin < 4,
        i = 0:varientNum-1;
    else
        i = floor(varientNum*rand(1,num));
    end
    fid = fopen(filename,'w');
    for i=i,
        fprintf(fid,'> Varient %d\n',i+1);
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