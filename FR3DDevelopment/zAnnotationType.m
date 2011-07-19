% zAnnotationType returns a number indicating Bulge, Helix, IL, HL, JL, etc.

% HL  -> 1
% H   -> 2
% IL  -> 3
% BL  -> 3   % treat it the same as an internal loop
% J   -> 5
% 3_p -> 6
% 5_p -> 7


function [a] = zAnnotationType(T)
  if strcmp(T(1:2),'HL'),
    a = 1;
  elseif strcmp(T(1:2),'IL'),
    a = 3;
  elseif strcmp(T(1:1),'J'),
    a = 5;
  elseif strcmp(T(1:1),'H'),
    a = 2;
  elseif strcmp(T(1:2),'BL'),
    a = 4;
  elseif strcmp(T(1:1),'5'),
    a = 6;
  elseif strcmp(T(1:1),'3'),
    a = 6;
  else
    T
    a = 7;
  end
