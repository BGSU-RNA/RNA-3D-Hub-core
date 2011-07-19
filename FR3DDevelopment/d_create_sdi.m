function [stacking_descrepancy_index stacking_descrepancy_matrix] = create_sdi(D)
stacking_descrepancy_index  = 0;
stacking_descrepancy_matrix = ones(4,4)/16;  % completely uninformative prior notion of replacement probabilities
y = D(:,1:2);                                % candidate bases
x = D(:,3);                                  % discrepancies
L = length(x);                               % length of discrepancy vector
w = zeros(L,1);                              % initialize unstandardized weighting vector of replacement possibilities
w_star = zeros(L,1);                         % initialize standardized weighting vector of replacement possibilities

 for i = 1:L
  w(i) = 1-2*x(i);                           % compute unstandardized weights on the interval (0,1), with higher weight
 end                                         %  given to smaller discrepancies.

 for i = 1:L
  w_star(i) = w(i)/sum(w);                   % standardize the weights to reflect emperical replacement probabilities.
 end

 for i = 1:L
  stacking_descrepancy_matrix(y(i,1),y(i,2)) = stacking_descrepancy_matrix(y(i,1),y(i,2)) + w_star(i); % weight the replacement probabilities
 end

 for i = 1:4
  for j = 1:4
   stacking_descrepancy_matrix(i,j) = stacking_descrepancy_matrix(i,j)/sum(sum(stacking_descrepancy_matrix)); % standardize the replacement probability matrix
  end
 end
 descrepancy_vector = sum(stacking_descrepancy_matrix');   
 stacking_descrepancy_index = std(descrepancy_vector);  % calculate the sdi like ICE.
end
