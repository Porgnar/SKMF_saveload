SKMF code capable of resuming, or continuing simulations.


FEATURES: 
-If you have a new type SKMF output file and a matching input file you can resume your simulation no matter how it quit
-If you have a new type SKMF output file and you know the radius it ran with, you can continue any simulation while changing the other parameters in the input file
-Freely change composition, parameters, timescale any time and resume your simulation, you can go from one equilibrium state to another one just by changing parameters
-Freely resume during the growth process, but keep in mind this way you may lose some accuracy on your average composition because a new set of random atoms get rolled for the yet to be deposited sites every time you resume


HOW TO USE:
-In the input.dat file change the Ns variable to the final step you want to run the simulation to.
-When starting the program if you give it no arguments (./SKMF) the simulation will completely overwrite any preexisting output.xyz file
-If you give the program a number argument (./SKMF 3149) the simulation will load and continue from the input step, running until the specified final step in input.dat


KNOWN BUGS IF YOU DON'T USE IT RIGHT
-If you resume from a step that doesnt exist, or is beyond the number of steps defined in input.dat the program will either segfault, or run but not do anything  