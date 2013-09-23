
some thoughts ...

- the CTEMgui should be capable of 
	reading output from any CTEM program and display it
	generating namelist input files for any CTEM program
	generating defect input files
	running the CTEM program in background and displaying runtime output

- input editor: shows all possible input fields, with those appropriate
  for a given program highlighted.  Program stores last used value for all
  fields and updates upon saving.

- program stores settings from previous run (folder, program, etc...)

- program logs each action

- three interfaces: nml input editor, defect editor, output display

- output display should have floating adjustable IMGWIN from Coyote's library

- for each CTEM program, we create a file that specifies the necessary inputs
  as well as the namelist name


