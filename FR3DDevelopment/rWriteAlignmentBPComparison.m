function [] = rWriteAlignmentBPComparison(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,sheet)

if nargin<6
  sheet=[];
end

[i j k] = find(File1.Edge(Indices1,Indices1));
BPAI=zeros(length(i),2);
ct=0;
for m = 1:length(i)
%     if (abs(k(m)) <= 20 || abs(k(m)) >=24) && (abs(k(m)) <= 120 || abs(k(m)) >= 124)
      if (abs(k(m)) < 14)
        ct=ct+1;
        BPAI(ct,1)=Indices1(i(m));
        BPAI(ct,2)=Indices1(j(m)); 
      end
end
BPAI=BPAI(1:ct,:);

[i j k] = find(File2.Edge(Indices2,Indices2));
BPBI=zeros(length(i),2);
ct=0;
for m = 1:length(i)
%     if (abs(k(m)) <= 20 || abs(k(m)) >=24) && (abs(k(m)) <= 120 || abs(k(m)) >= 124)
    if (abs(k(m)) < 14)
        ct=ct+1;
        BPBI(ct,1)=Indices2(i(m));
        BPBI(ct,2)=Indices2(j(m));
    end
end
BPBI=BPBI(1:ct,:);

BPList1=zeros(length(BPAI)+length(BPBI),2);
BPList2=zeros(length(BPAI)+length(BPBI),2);
ct=0;
for i=1:length(BPAI)
   BPList1(i,1)=BPAI(i,1);
   BPList1(i,2)=BPAI(i,2);
   ct=ct+1;
   p=find(AlignedIndices1==BPAI(i,1));
   r=find(AlignedIndices1==BPAI(i,2));
   if isempty(p)
      BPList2(i,1)=0;
      if isempty(r)
         BPList2(i,2)=0;
      else
         BPList2(i,2)=AlignedIndices2(r);
      end
   else
      if isempty(r)
         BPList2(i,1)=AlignedIndices2(p);
         BPList2(i,2)=0;
      else
         BPList2(i,1)=AlignedIndices2(p);
         BPList2(i,2)=AlignedIndices2(r); 
         q=find(BPBI(:,1)==AlignedIndices2(p));
         s=find(BPBI(:,2)==AlignedIndices2(r));
         t=intersect(q,s);
         if ~isempty(t)
            BPBI(t,:)=[];
         end
      end
   end
end
ct=length(BPAI);
for i=1:length(BPBI)
   ct=ct+1;
   BPList2(ct,1)=BPBI(i,1);
   BPList2(ct,2)=BPBI(i,2);
   p=find(AlignedIndices2==BPBI(i,1));
   r=find(AlignedIndices2==BPBI(i,2));
   if isempty(p)
      BPList1(ct,1)=0;
      if isempty(r)
         BPList1(ct,2)=0;
      else
         BPList1(ct,2)=AlignedIndices1(r);
      end
   else
      if isempty(r)
         BPList1(ct,1)=AlignedIndices1(p);
         BPList1(ct,2)=0;
      else
         BPList1(ct,1)=AlignedIndices1(p);
         BPList1(ct,2)=AlignedIndices1(r); 
         q=find(BPAI(:,1)==AlignedIndices1(p));
         s=find(BPAI(:,2)==AlignedIndices1(r));
         t=intersect(q,s);
         if ~isempty(t)
            BPAI(t,:)=[];
         end
      end
   end
end
for i = 1:length(AlignedIndices1)
   if ~any(BPList1(:,1)==AlignedIndices1(i)) || ~any(BPList2(:,1)==AlignedIndices2(i))
      ct=ct+1;
      BPList1(ct,1)=AlignedIndices1(i);
      BPList1(ct,2)=0;
      BPList2(ct,1)=AlignedIndices2(i);
      BPList2(ct,2)=0;
   end
end
for i = 1:length(Indices1)
   if ~any(BPList1(:,1)==Indices1(i))
      ct=ct+1;
      BPList1(ct,1)=Indices1(i);
      BPList1(ct,2)=0;
      BPList2(ct,1)=0;
      BPList2(ct,2)=0;
   end
end
for i = 1:length(Indices2)
   if ~any(BPList2(:,1)==Indices2(i))
      ct=ct+1;
      BPList1(ct,1)=0;
      BPList1(ct,2)=0;
      BPList2(ct,1)=Indices2(i);
      BPList2(ct,2)=0;
   end
end
BPList1=BPList1(1:ct,:);
BPList2=BPList2(1:ct,:);
[V,I]=sort(BPList1(:,1),'ascend');
BPList1=BPList1(I,:);
BPList2=BPList2(I,:);

num2move=length(find(BPList1(:,1)==0));
N=length(BPList1(:,1));
for i=1:num2move
   tmpBP1=BPList1(1,:);
   tmpBP2=BPList2(1,:);
   P=find(BPList2(num2move+1:N,1)>tmpBP2(1,1),1,'first');
   if isempty(P)
      BPList1(1:N-1,:)=BPList1(2:N,:);
      BPList1(N,:)=tmpBP1;
      BPList2(1:N-1,:)=BPList2(2:N,:);
      BPList2(N,:)=tmpBP2;
   else
      BPList1(1:num2move+P-2,:)=BPList1(2:num2move+P-1,:);
      BPList1(num2move+P-1,:)=tmpBP1;
      BPList2(1:num2move+P-2,:)=BPList2(2:num2move+P-1,:);
      BPList2(num2move+P-1,:)=tmpBP2;
   end
   num2move=num2move-1;
end

FinalListing=cell(length(BPList1(:,1)),6);
for i=1:length(FinalListing)
     if BPList1(i,1)==0
         a='---';
     else
         a = File1.NT(BPList1(i,1)).Number;
     end
     if BPList1(i,1)==0 || BPList1(i,2)==0
        b=' ';
     else
         b = zEdgeText(File1.Edge(BPList1(i,1),BPList1(i,2)));
     end
     if BPList1(i,2)==0
         c='---';
     else
         c = File1.NT(BPList1(i,2)).Number;
     end
     if BPList2(i,1)==0
         d='---';
     else
         d = File2.NT(BPList2(i,1)).Number;
     end
     if BPList2(i,1)==0 || BPList2(i,2)==0
        e=' ';
     else
        e = zEdgeText(File2.Edge(BPList2(i,1),BPList2(i,2)));
     end
     if BPList2(i,2)==0
         f='---';
     else
         f = File2.NT(BPList2(i,2)).Number;
     end
     
     FinalListing{i,1}=a;
     FinalListing{i,2}=b;
     FinalListing{i,3}=c;
     FinalListing{i,4}=d;
     FinalListing{i,5}=e;
     FinalListing{i,6}=f;
end

ExcelName=[File1.Filename '_' File2.Filename '_Alignment']; 

if isempty(sheet)
   xlswrite(ExcelName,FinalListing)
else
   xlswrite(ExcelName,FinalListing,sheet)
end
Excel = actxserver('Excel.Application');
% eWorkbook = Excel.Workbooks.Add;
Excel.Workbooks.Open(ExcelName);

for j=[2 5]
   for i=1:length(FinalListing)
      if j==2
         Range = Excel.Range(['B' num2str(i)]);
      else
         Range = Excel.Range(['E' num2str(i)]);
      end
      switch lower(strtrim(FinalListing{i,j}))
         case 'cww'
            Range.Interior.ColorIndex = 4;
         case {'tsh','ths'}
            Range.Interior.ColorIndex = 38;
         case {'cws','csw'}
            Range.Interior.ColorIndex = 33;
         case {'tws','tsw'}
            Range.Interior.ColorIndex = 41;     
         case 'tss'
            Range.Interior.ColorIndex = 13;
         case {'thw','twh'}
            Range.Interior.ColorIndex = 45; 
         case 'css'
            Range.Interior.ColorIndex = 54;
         case {'chs','csh'}
            Range.Interior.ColorIndex = 7; 
         case {'chw','cwh'}
            Range.Interior.ColorIndex = 40;
         case 'thh'
            Range.Interior.ColorIndex = 12;
         case 'tww'
            Range.Interior.ColorIndex = 6;   
         case 'chh'
            Range.Interior.ColorIndex = 9;
      end
   end
end
% Excel.Workbooks.Saved = 1;
% Excel.Workbooks.saveas(ExcelName);
% Excel.Workbooks
% eWorkbook
% eWorkbook.SaveAs([ExcelName 'v2']);
% eWorkbook.Saved = 1;
Excel.Workbooks.Close;
Excel.Quit;
Excel.delete;
