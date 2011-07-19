

% load a structure

% F = zAddNTData('1s72');
Indices = zIndexLookup(F,'75:81_9,101:106_9');
% Letters = GGAGUAC*GGAAAC

Variant{1}     = 'CAGUAGAAC';
Variant{16}    = 'CAGUAGACA';
Variant{41}    = 'CAGUAGAAG';
Variant{42}    = 'UAGUUGAAA';
Variant{84}    = 'CAGUUGACA';
Variant{118}   = 'GAGUAAAAA';
Variant{310}   = 'CACUAGAAC';
Variant{521}   = 'CAACAGACA';
Variant{8530}  = 'GGGUAAAAA';
Variant{96454} = 'AGAAUAAAU';


for v = 1:length(Variant),
  if ~isempty(Variant{v}),
    V = Variant{v};
    Letters = ['G' V(1:5) 'CG' V(6:9) 'C'];
    Motif = xSubstituteNucleotides(F,Indices,Letters,1);
    pause
    zWritePDB(Motif,['Sarcin_' Letters '_Variant_' num2str(v) '.pdb']);
  end
end
