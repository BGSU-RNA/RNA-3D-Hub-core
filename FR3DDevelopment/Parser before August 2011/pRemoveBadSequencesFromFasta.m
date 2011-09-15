
%Sequence = pReadFASTA('5S_Rfam_Archaea_seed_Jesse_2_20_05.fasta',1,Inf);
Sequence = pReadFASTA('2004_16sB_UNIQUE.fasta',1,Inf);

 % --------------- exclude sequences with letters other than ACGU.-

  for k=1:length(Sequence),
    NumOther(k) = sum(1-ismember(Sequence(k).X,'ACGU.-'));
    Length(k)   = length(Sequence(k).X);
  end
  figure(1)
  clf
  hist(NumOther,30)

  figure(2)
  clf
  hist(Length,30)

  i = find((NumOther < 76) .* (Length > 1400) .* (Length < 1800)); 
  j = find((NumOther >= 76) +  (Length <= 1400) +  (Length >= 1800)); 

  GoodSequence = Sequence(i);                       % remove those w/ difft letters

  pWriteFasta(GoodSequence,'2004_16sB_UNIQUE_Less_than_76_non_ACGU.fasta');

  BadSequence = Sequence(j);                       % remove those w/ difft letters

  pWriteFasta(BadSequence,'2004_16sB_UNIQUE_Rejected.fasta');


