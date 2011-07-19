% Produce a text output in N3 format of the FR3D-classified interactions in the given RNA molecule

function [T] = zTextRNAOAssertions(File)

r = 1;

T{r} = 










@prefix rna: <http://purl.obolibrary.org/obo/>.
@prefix ro: <http://purl.obolibrary.org/obo/>.
@prefix : <http://example.org>.
:m1 rdfs:label "".
:m1 :pdb_code "1A50"
:m1 rdf:type rnao:RNAO_0000168.                                   # m1 is a molecule or something
:m1 rnao:has_proper_part :nt1.                                        # m1 contains something called nt1
:nt1 rdf:type rnano:RNAO_0000102.                                 # nt1 is a cytosine residue
:m1 rnao:has_proper_part :nt2.                                        # m1 contains something called nt1
:nt2 rdf:type rnao:RNAO_0000105.                                   # nt2 is a guanine residue
:nt1 rnao:three_prime_to_five_prime_to :nt2
:nt1 rnao:pairs_with_CSW :nt2
