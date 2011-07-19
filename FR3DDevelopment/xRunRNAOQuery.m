
Query = xReadRNAOQuery('RNAOQuery.txt',1);
Query.SearchFiles = {'1s72'};
Query.SearchFiles = {'Nonredundant_4A_2010-05-19_list'};
Filenames = Query.SearchFiles;
Query = xConstructQuery(Query);
UsingLibrary = 1;
xFR3DSearch
xListCandidates(Search)
xDisplayCandidates(File,Search)
