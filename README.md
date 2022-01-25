# Generalized-Replica-Exchange
A python wrapper for setting up, running, and analyzing gREM simulations with LAMMPS.
This code is made up of multiple parts that can be used to run and analyze simulations using the generalized replica exchange method.

This code is made up of multiple parts that can be used to run and analyze simulations using the generalized replica exchange method.

The main part of the code is “run_gREM.py,” which uses some classes that are defined in other files (“sortdumps.py”, and “walker.py”)


In this readme file, I will briefly describe how to setup, run, and analyze a GREM simulation using these codes.

Things that are assumed:
1) You want lambdas that are evenly spaced
2) You know the minimum, maximum, and spacing of those lambdas that you want (if you don’t, read the procedure from Kim et al that describes the GREM method).



##Setup Stage

To setup a GREM simulation, you need to use the setup option of the GREM code.

You can do this by running:

run_gREM.py -setup 1 -start 300 -stop 1000 -nbins 71

Here, the setup option overrides all later parts of the code and ONLY sets up the directory structure and some basic files. 

The start option in this case is the lowest lambda value. 

The stop option in this case is the highest lambda value

The nbins option is the number of windows/walkers that you want to run.

Lambdas (and their spacing) is selected by evenly distributing nbins between start and stop.



##Walkdown Stage

The next step is to run a “walkdown” where you run enthalpies in order from lowest lambda to highest lambda and launch each successive lambda from the previous lambda. This makes sure that you are equilibrating from the lowest state.

run_gREM.py -start 0 -stop 71 -eta -2.5 -Ho -1800 -nprocs 16 -inp walkdown.in -lambdafile grem.include -walkdown 1

Here start and stop now serve the purpose of identifying the range of windows to use for the walkdown.

Eta refers to the slope of the function determining temperature

Ho is the reference enthalpy

Nprocs is the number of processors to use

Inp is the lammps input file for the walkdown

Lambdafile is the file generated in the setup phase which includes the lammps variable defs for lambda and replica

Walkdown 1 turns on the walkdown step.

Once you run this, all you have to do is submit the submission script as qsub walkdown.sh - if all is well every trajectory will run successfully. You are now ready to start actual gREM.


##Run GREM using LAMMPS

From this point, the next step is to run actual gREM using lammps with existing scripts.


##Analyze GREM Runs

To do this, all you have to do is run without setting the walkdown or setup options. It is suggested you make a directory elsewhere and use the workdir to point to the right path as this step generates tons of files.

Start/Stop define now the range of runs (if there are more than one) to use. The first time you run you will want to use the restart 0 [default] option, which just reads everything into python pickle files for fast analysis

Nbins determines the number of histogram bins for use with ST-WHAM

When you rerun the analysis code after running once with restart 0, and flip it to restart 1, it then will read the restart files, do ST-WHAM, and any other final analysis as needed.

##Example:

run_gREM.py -start 1 -stop 6 -eta -0.1 -Ho -19600 -workdir .. -nbins 50 -nb 5 -restart 0
run_gREM.py -start 1 -stop 6 -eta -0.1 -Ho -19600 -workdir .. -nbins 50 -nb 5 -restart 1

The first does initial setup, and the second does final analysis.

dump_analyze.py -dumpid 0 -dumpbase dump -nb 5 -start 1 -stop 6 -workdir .. -nreps 86 -option 1


Changelog


January 11, 2022 Added an option to dump a pickle file with the block data for H. Shape = (#replicas, #blocks)
