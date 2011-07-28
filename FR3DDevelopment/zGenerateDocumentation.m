
disp('Make BGSURNA the working directory');

folder = 'MutationRates';
docfolder = [folder 'Doc'];
if exist(docfolder) == 7,
  rmdir(docfolder,'s');
end
m2html('mfiles',folder, 'htmldir',docfolder, 'graph','on');

folder = 'FR3DSource';
docfolder = [folder 'Doc'];
if exist(docfolder) == 7,
  rmdir(docfolder,'s');
end
m2html('mfiles',folder, 'htmldir',docfolder, 'graph','on');

folder = 'FR3DDevelopment';
docfolder = [folder 'Doc'];
if exist(docfolder) == 7,
  rmdir(docfolder,'s');
end
m2html('mfiles',folder, 'htmldir',docfolder, 'graph','on');


