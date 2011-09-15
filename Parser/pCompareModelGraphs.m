function pCompareModelGraphs(SeqFile,ModFile1,ModFile2)
    NumSeq=1000;
    PD1 = JAR3DMatlab.Align2(pwd,SeqFile,ModFile1,NumSeq,0,15);
    PD2 = JAR3DMatlab.Align2(pwd,SeqFile,ModFile2,NumSeq,0,15);
    probsM1 = double(PD1.probsM(1:size(PD1.probsM),1:size(PD1.probsM(1))));
    probsM2 = double(PD2.probsM(1:size(PD2.probsM),1:size(PD2.probsM(1))));
    y1=mean(probsM1,1);
    y1(1:length(y1)-1)=y1(1:length(y1)-1)-y1(2:length(y1));
    y1
    y2=mean(probsM2,1);
    y2(1:length(y2)-1)=y2(1:length(y2)-1)-y2(2:length(y2));
    y2
    x1=(0:(length(y1)-1))/(length(y1)-1);
    x2=(0:(length(y2)-1))/(length(y2)-1);
    plot(x1,y1,x2,y2);
end
