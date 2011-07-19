
function [Text] = zBasepairCodeText(b)

switch b
  case  2, Text = 'cWW';
  case  3, Text = 'tWW';
  case  4, Text = 'cWH';
  case  5, Text = 'tWH';
  case  6, Text = 'cWS';
  case  7, Text = 'tWS';
  case  8, Text = 'cHH';
  case  9, Text = 'tHH';
  case 10, Text = 'cHS';
  case 11, Text = 'tHS';
  case 12, Text = 'cSs';
  case 13, Text = 'tSs';
  case 14, Text = 'bif';
  case 15, Text = 'cHW';
  case 16, Text = 'tHW';
  case 17, Text = 'cSW';
  case 18, Text = 'tSW';
  case 19, Text = 'cHH';
  case 20, Text = 'tHH';
  case 21, Text = 'cSH';
  case 22, Text = 'tSH';
  case 23, Text = 'csS';
  case 24, Text = 'tsS';
  case 25, Text = 'bif';
end
