
Each resource file (something.res) has as its first line the 
number of namelists that are used by the program (typically 1)
followed by the name of each namelist block, the number of 
namelist entries and a numerical code mapping the entry to
a master list of entries.  So, the template for these files
is as follows

2				% # of namelists
namelistname1			% name of namelist
3				% number of entries in this namelist
200::name			% entry ID in masterlist
600::name
700::name
namelistname2			% etc ...
200::name
203::name
701::name


The file DEFAULTS.txt shows for each numerical index the corresponding 
variable name and a default value in the format '###::name::value', one
line per variable.


The file SNAPSHOT.txt is initially identical to the DEFAULTS.txt file,
but will change as it keeps the last user entered value for each entry.



