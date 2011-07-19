% makeRandomSequencesFasta(filename,dim)

function [Sequence] = makeRandomSequencesDLFasta(filename,dim)

letter = ['A','C','G','U'];
cWW = ['AU','CG','GC','UA','GU','UG'];
fid = fopen(filename,'w');
[strands,N]=size(dim);
Sequences{N} = '';                      % allocate space for all sequences
if strands == 2,
    for i = 1:N,                                % N sequences
        LL = max(dim(1,i),2);
        RL = max(dim(2,i),2);
        fprintf(fid,'> Variant %d\n',i);
        openp = (ceil(6*rand(1,1)));
        closep = (ceil(6*rand(1,1)));
        fprintf(fid,'%s',cWW(2*openp-1));
        s = cWW(2*openp-1);
        m = ceil(4*rand(1,LL-2));               % uniform distribution ACGU
        fprintf(fid,'%s',letter(m));
        s = [s letter(m)];
        fprintf(fid,'%s*%s',cWW(2*closep-1),cWW(2*closep));
        s = [s cWW(2*closep-1) '*' cWW(2*closep)];
        m = ceil(4*rand(1,RL-2));               % uniform distribution ACGU
        fprintf(fid,'%s',letter(m));
        s = [s letter(m)];
        fprintf(fid,'%s\n',cWW(2*openp));
        s = [s cWW(2*openp)];
        Sequence{i} = s;
    end
else
    for i = 1:N,
        fprintf(fid,'> Variant %d\n',i);
        openp = (ceil(6*rand(1,1)));
        fprintf(fid,'%s',cWW(2*openp-1));
        m = ceil(4*rand(1,dim(i)-2));           % uniform distribution ACGU
        fprintf(fid,'%s',letter(m));
        s = [cWW(2*openp-1) letter(m) cWW(2*openp)];
        fprintf(fid,'%s\n',cWW(2*openp));
        Sequence{i} = s;
    end
end
fclose(fid);
end
