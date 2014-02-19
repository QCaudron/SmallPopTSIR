Utilities
=========


Files
-----

*uk_generator.py*

Takes the data from /raw_data/uk and generates /data/uk_all containing a .csv for each town in the UK. Ensures 
interpolations and normalisations are the same everywhere.

*uk_stats.py*

Runs simulations for each .csv in /data/uk_all, and collects data on predicted vs real sizes and durations, as 
well as S0 vs final size, into one .csv. Output files are in /data/uk_all/generated/. 

*notebook.ipynb*

Like batch.py, but in IPYNB format for development.

*data_protocol.jpg*

What we did with the raw data to get useable data for simulation.
