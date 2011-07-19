%----------------------------------------------------------------------------

function [s] = zGetAtomStorageLocation(name,c);

switch name,
  case 'C1*',     s = 'Sugar(1,:)';
  case 'C2',      if (c==1) | (c==3), s = 'Fit(9,:)'; else s = 'Fit(2,:)'; end
  case 'C2*',     s = 'Sugar(2,:)';
  case 'C3*',     s = 'Sugar(4,:)';
  case 'C4',      if (c==1) | (c==3), s = 'Fit(2,:)'; else s = 'Fit(5,:)'; end
  case 'C4*',     s = 'Sugar(6,:)';
  case 'C5',      s = 'Fit(8,:)';
  case 'C5*',     s = 'Sugar(8,:)';
  case 'C6',      if (c==1) | (c==3), s = 'Fit(5,:)'; else s = 'Fit(7,:)'; end
  case 'C8',      s = 'Fit(7,:)';
  case 'H1',      if     (c==2), s = 'Fit(9,:)'; 
                  elseif (c==3), s = 'Fit(12,:)'; 
                  else s = 'Fit(10,:)'; end
  % A has no H1 atom
  case 'H2',      s = 'Fit(11,:)';
  case 'H3',      s = 'Fit(11,:)';
  case 'H5',      if     (c==2), s = 'Fit(11,:)'; 
                  else s = 'Fit(9,:)'; end
  case 'H6',      if     (c==2), s = 'Fit(10,:)'; 
                  else s = 'Fit(12,:)'; end
  case 'H8',      if     (c==1), s = 'Fit(12,:)'; 
                  else s = 'Fit(13,:)'; end
  case 'H9',      if     (c==1), s = 'Fit(13,:)'; 
                  else s = 'Fit(14,:)'; end
  case '1H1',     s = 'Fit(14,:)';
  case '1H6',     s = 'Fit(14,:)';
  case '2H6',     s = 'Fit(15,:)';
  case '1H4',     s = 'Fit(12,:)';
  case '2H4',     s = 'Fit(13,:)';
  case '1H2',     s = 'Fit(15,:)';
  case '2H2',     s = 'Fit(16,:)';
  case 'N1',      if (c==1) | (c==3), s = 'Fit(4,:)'; else s = 'Fit(1,:)'; end
  case 'N2',      s = 'Fit(11,:)';
  case 'N3',      if (c==1) | (c==3), s = 'Fit(3,:)'; else s = 'Fit(4,:)'; end
  case 'N4',      s = 'Fit(6,:)';
  case 'N6',      s = 'Fit(6,:)';
  case 'N7',      s = 'Fit(10,:)';
  case 'N9',      s = 'Fit(1,:)';
  case 'O1P',     s = 'Sugar(11,:)';
  case 'O2',      s = 'Fit(3,:)';
  case 'O2*',     s = 'Sugar(3,:)';
  case 'O2P',     s = 'Sugar(12,:)';
  case 'O3*',     s = 'Sugar(5,:)';
  case 'O4',      s = 'Fit(6,:)';
  case 'O4*',     s = 'Sugar(7,:)';
  case 'O5*',     s = 'Sugar(9,:)';
  case 'O6',      s = 'Fit(6,:)';
  case 'P',       s = 'Sugar(10,:)';
  otherwise,      s = '????????'; fprintf('%s\n',name);
end
