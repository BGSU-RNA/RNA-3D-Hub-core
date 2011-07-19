
function [Text] = zInteractionCodeText(b)

Text = '   ';

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

  case 30, Text = 's35';
  case 31, Text = 's53';
  case 32, Text = 's33';
  case 32, Text = 's33';
  case 33, Text = 's55';
  case 33, Text = 's55';

  case 40, Text = 'ncWW';
  case 41, Text = 'ntWW';
  case 42, Text = 'ncWH';
  case 43, Text = 'ntWH';
  case 44, Text = 'ncWS';
  case 45, Text = 'ntWS';
  case 46, Text = 'ncHH';
  case 47, Text = 'ntHH';
  case 48, Text = 'ncHS';
  case 49, Text = 'ntHS';
  case 50, Text = 'ncSS';
  case 51, Text = 'ntSS';
  case 52, Text = 'nbif';
  case 53, Text = 'ncHW';
  case 54, Text = 'ntHW';
  case 55, Text = 'ncSW';
  case 56, Text = 'ntSW';
  case 57, Text = 'ncHH';
  case 58, Text = 'ntHH';
  case 59, Text = 'ncSH';
  case 60, Text = 'ntSH';
  case 61, Text = 'ncSS';
  case 62, Text = 'ntSS';
  case 63, Text = 'nbif';
  case 70, Text = 'ns35';
  case 71, Text = 'ns53';
  case 72, Text = 'ns33';
  case 73, Text = 'ns55';
end
