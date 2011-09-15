function pOneModelTest(SeqFile,ModFile)
    NumSeq=1000;
    PD1 = JAR3DMatlab.Align2(pwd,SeqFile,ModFile,NumSeq,0,15);
    probsM1 = double(PD1.probsM(1:size(PD1.probsM),1:size(PD1.probsM(1))));
    if(min(size(probsM1)>1))
        y1=mean(probsM1,1);
    else
        y1=probsM1;
    end
    probsM1
    % y1(1:length(y1)-1)=y1(1:length(y1)-1)-y1(2:length(y1));
    y1
    x1=(0:(length(y1)-1))/(length(y1)-1);
    plot(x1,y1);
end