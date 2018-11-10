# Ripple detection with Trodes analysis and evaluation code
<pre>Update coming soon (11-10-18)</pre>

Note: `nelpy` is widely used during this analysis. However, the version of `nelpy` used is a modified one. See links below:<br><br>
Nelpy<br>
[Nelpy used for analysis](https://github.com/shayokdutta/nelpy_modified)<br>
[Nelpy github](https://github.com/nelpy)
<br><br>
Trodes<br>
[Trodes + RippleDetection module software](https://bitbucket.org/mkarlsso/trodes/branch/rippleDetectionBeagleBoneStimModule)<br>
[How to use Trodes + RippleDetection module](https://docs.google.com/document/d/1cZG8eLavlUdqzHJkOljiraBl5P1K1ZV-zJhWvluX0Tg/edit?usp=sharing)
<br><br>
Paper<br>
[bioRxiv preprint](https://www.biorxiv.org/content/early/2018/04/11/298661)
<br><br>
Posters and other links can be found on [sd54.web.rice.edu](http://sd54.web.rice.edu/research/).

## File Descriptions
This is being updated. This is not completed yet<br><br>

### General File(s)
1. `DataAnalysisScripts/ripple_filtering.py` basic ripple detection functions<br>

### Single Channel Detection Analysis
2. `DataAnalysisNotebooks/SingleChannel.ipynb` has canonical single channel ripple detection

3. `DataAnalysisScripts/singleChannelSimulatedDetections/singleChanAnalysis_singleChanDefn.cpp` varies detection threshold params while simulating single channel detections<br>

4. `DataAnalysisScripts/singleChannelSimulatedDetections/singleChanAnalysis_twoChanDefn.cpp` varies detection threshold params while simulating single channel detections on two channel cannonical ripple definition<br>

### Multi-channel Detection Analysis
5. `DataAnalysisNotebooks/MultichannelAnalysis.ipynb` has canonical detections using two channels<br>

6. `DataAnalysisScripts/voting/voting2of2.cpp` varies detection threshold params while simulating ripple detections across 2 channels that have to agree within 15 ms<br>

7. `DataAnalysisScripts/voting/voting2of3.cpp` varies detection threshold params while simulating ripple detections across 3 channels in which 2 of which agree that ripples are detected within 15 ms of each other <br>

### Synthetic Ripple Detection Analysis
8. `DataAnalysisNotebooks/SyntheticRippleAnalysis.ipynb` has general pipeline for generating synthetic ripples and explores parameter space across varying synthetic ripples<br>

9. `DataAnalysisScripts/simulatedRippleDetectionsSingleChan/main.cpp` performs simulated ripple detections on synthetic data using a single channel with varying ripple amplitudes and 

## How to modify and run things
This is not complete. Unfortunately, a lot if this requires you to figure it out on your own. I shall make it more complete over time but if you want to explore before it has been updated feel free to reach out to me.<br><br>
In general, all analysis begins with the filtering/thresholding functions located in the file listed after 1 in the File Descriptions section. In order to run this with your own data or own filter coefficients and run the evaluations that we have run in the paper, follow the steps below. <br><br>

1. Modify the filter coefficients within the `ripple_filtering.py` file.<br>

2. Grab a `.rec` file from `Trodes` or your own data and load it up into `Python` with the same filtering that is done during the realtime detections or expected to happen during the realtime detections. E.g., the version of `Trodes` used in our RippleDetector has a 400 Hz lowpass filter for LFP so data is first lowpass filtered at 400 Hz with modified extraction functions provided within the bitbucket repository for the RippleDetector. <br>

3. Perform ripple band filtering, power signal calculations and threshold crossing detections in order to determine canonical ripples. The code for these are provided within `ripple_filtering.py`. Examples of how to use this are within `SingleChannel.ipynb` and `MultichannelAnalysis.ipynb`. Save the canonical detection start and end bounds in separate files. See files 3, 4, 6 and 7 from the File Descriptions section to examine how they are being processed. <br>

4. Perform ripple band filtering and power estimations for how it will be in the realtime detections. Examples of these are shown in the analyisis notebooks. Save the final power estimation signal (i.e., the signal after the ripple band filtering, absolute value and lowpass filtering).<br>

5. Modify realtime detection estimation files (voting cpp files and/or singleChanAnalysis cpp files) to have desired file names that are defined in the macros at the top. Run the realtime detection code once it is setup appropriately!<br>

6. To visualize and plot results, modify the files within the `figureCode_final` directory. In particular, `Figure4_SingleChannelAnalysis.ipynb` and `Figure5_MultiChannelAnalysis.ipynb`. How the outputs should look in order to input them into these files are provided within the `paperData/figureData` file directory.
