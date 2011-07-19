% zCountTriples is used to load a SearchSaveFile from a search for triples, then count how many unique triples appear there, and how many 2-basepair and 3-basepair triples as well

% load 2010-07-11_00_08_21-Coplanar_triples
% load 2010-07-20_10_26_34-Triples_in_23S
% load 2010-07-20_10_30_13-Triples_in_16S
% load 2010-07-20_10_34_59-pair_cp_ncp_not_stack_triples
% load 2010-07-20_10_54_01-pair_cp_ncp_not_stack_one_not_pair_triples
% load 2010-07-20_10_56_15-cp_ncp_not_stack_no_pair_triples
% load 2010-07-20_11_45_16-cp_ncp_not_stack_triples_excl_1QCU
load 2010-07-20_10_32_36-cp_ncp_not_stack_triples

C = Search.Candidates;
C(:,1:3) = sort(C(:,1:3),2);
  
[U,n] = zUniqueRows(C);
n = n';
[V,m] = zUniqueRows(n);

for v = 1:length(V),
  fprintf('There were %4d sets listed %2d times.\n', m(v), V(v)); 
end

fprintf('Based on this, there are %4d distinct triples.\n',length(n));
