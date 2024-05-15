README.md
### This is a living document of changes and updates to code and datasets ### 

Repo Folders

code: 

birdnet_data_process.R 
- This script is for converting birdnet output from a bunch of .txt files to .csv format
- Once in .csv format you can then

       - Add alpha codes, scientific names, format columns for date and time
       - filter down to specific sites
       - filter down to specific times
       - filter by confidence score (expanding this out soon)
       - check for gaps in data by date
       - look at number of unique species detected and total counts

data: 

data_2023_07_output.zip 
- This is a folder of unprocessed .txt files produced by birdnet from all recordings in July of 2023 
- be careful it is a lot of data
  
acoustic_dat_07.csv 
- this is basically all the files from data_2023_07_output.zip strung together and converted to a .csv

r3_07_23.csv 
- All diurnal observations from site "R3" with some other small modifications

species.list.csv 
- list of common names, alpha codes, and scientific names of birds 

Updates

updates for week of 05/13-05/17/2024
- Getting things setup. pardon the mess
