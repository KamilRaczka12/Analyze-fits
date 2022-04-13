# Analyze-fits
Program to analyze fits through Astrometry, Source Extractor and find it in NASA horizons
You need: Astrometry, Source Extractor, Python3 with packages specified in the provided 'requirements.txt' file
Make sure to check absolute paths in your sextractor's 'default.sex' file and astrometry's 'astrometry.cfg' file 
In 'default.sex' change 'default.param' to provided 'parameters.param' file with absolute path to that file
For FitsToFiles.py file:
you need to know your exact headers for RA, DEC and DATE from your fits files
example of how RA looks like: '06:43:22.1' or '09 01 47.7'
program accounts for both versions
you need to add path to folder with fits files and folder for sextractor config
