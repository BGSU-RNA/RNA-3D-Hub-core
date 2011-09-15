% pModelPDF calculates the percentage of scores in the empirical
% distribution for the model's  scores are lower than the scores given in
% the vector of test scores

% inputs  : test_scores; a vector of test scores from the given model, and
%           model; the model number for the model of interest, given as int
% outputs : p_vec; a vector of the same length as test scores giving the 
%           percentage of the scores in the empirical distribution less
%           than than each score in test_scores
%           mean_p; the percentage of the scores in the empirical distribution less
%           than than the mean score in test_scores

function [p_vec,mean_p] = pModelPDF(test_scores, model)
    if model < 10
        ms = ['00' int2str(model)];
    elseif model < 100
        ms = ['0' int2str(model)];
    else
        ms = int2str(model);
    end
    distFN = ['Models' filesep 'Emp. Distributions' filesep 'IL_' ms '.mat'];
    load(distFN,'Scores','-mat');
    mean_test_score = mean(test_scores);
    mean_p = mean(Scores <= mean_test_score);
    n = length(test_scores);
    p_vec = zeros(0,n);
    for i = 1:n
        p_vec(i) = mean(Scores <= test_scores(i));
    end
end