# TimeSeriesProcess-CZI
Automated processing of CZI data obtained on Zeiss confocal microscope; for time series experiments made up of line scans.

Designed to process time series data collected on a Zeiss confocal microscope from microfluidic flow experiments.


Data collection: A line scan is selected at one point along the microfluidics channel; the line scan is repeated to create a time series experiment. Time series experiments are repeated for a control (in which the fluorescence intensity profile should be flat), followed by 2 experiments in which the fluorescence intensity profile is expected to be perturbed. The experiments are obtained at different flow rates, to lead to different residence times at the same point along the microfluidics channel.

Data storage: Each flow rate has it’s own directory, inside which the 3 time series experiments (control, and two others) are saved as CZI files. The control always appears at the top of the direction, alphabetically, by appropriate file naming (important for processing). 

Using script: Run Cell 0 to set up processing for all flow rates (defining lists: dsets and FR is important for comparing across flow rates); only run once at the beginning. Import and check the averaged data across each time series in a single directory (flow rate) by running Cell 1, and normalizing with Cell 2. If the data is acceptable, use Cell 3 to store the data and flow rate (in dsets and FR, respectively). Repeat Cells 1-3 for each directory. Once all of the directories of interest have been imported, checked, normalized and stored (using Cells 1-3), run Cell 4 to plot a comparison of the results across directories, which in this case comprises of results at different flow rates. 
