SmallPopTSIR
============

Predictability of measles epidemics in small populations, using TSIR model fitting and simulation. This 
repository goes with the paper :


Caudron Q, Mahmud A S, Metcalf C J E, Gottfredsson M, Viboud C, Cliff A D, Grenfell B T, 
*Predictability in a highly stochastic system : final size of measles epidemics in small populations*. Journal of the Royal Society Interface **12** (102), 2014.




Files and Directories
---------------------

**batch.py** *directory_name number_of_simulations* 

Main Python script for processing an entire directory of data. Arguments : `directory_name number_of_sims` where `directory_name` is a directory in /data, and `number_of_sims` is an integer number of simulations to run for each time series. Call `python batch.py iceland 1000` for one thousand simulations for each time-series in Iceland, for example.



**data/**

Processed data in a form useable by all scripts, and hopefully ssm at some point.



**poster/**

Poster for the Ecology and Evolution of Infectious Diseases 2014 conference.



**raw_data/**

All raw data collected from stats agencies ( or Bryan Grenfell, or others ) go here for permanent storage.


**talk/**

Talk for Tim Coulson's group at Oxford's Department of Zoology.


**utilities/**

Random scripts used in the development of the code or in preprocessing of data.
