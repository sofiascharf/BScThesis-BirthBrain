# Bachelors-BirthingBrain
A GitHub Repository containing all final code for my bachelor's thesis 
The repository is divided into the following folders: 

## Raw-Data 
A folder containing raw data used for analysis 

### Audio Data
.wav files recorded with USV microphones, through AVISOFT software.

### Video Data 
.mp4 files recorded with USB cameras, through EthoVisionXT.

### Fiber Photometry Data
.??? files recorded with optic fibres, through DoricStudio



## On-line-Detection-USV
A script used for on-line detection parturient behavior through USV detections. 
The script runs infinitely ones you press start, and you have to manually start the USV recordings right before.
The script will be divided into different sections for collecting, processing and predicting the data. 
This script is created and updated throughout based on findings in another analysis section:



## Behavior-Analysis

## USV Analysis 
This contains a folder with the Detection Files created through DeepSqueak, and the code for creating them.
It further contains a script for concatenizing separate Detection Files from one session into a joint .csv file, and that said .csv file.
Finally, it contains two notebooks. 

### Notebooks
One notebook for *Data Visualization*: This is used to explore the data and visually inspect for differences that can mark pre-parturient behavior. This involves converting the .csv file into time series data. Furthermore, I might compute metrics to support a hypothesis that mice exhibit dissociation behavior. 

One notebook for *Classification*: This is used to create the classifier used for on-line detection. Different types of classifiers are compared, and different options such as choosing the right baseline recordings, and optimal temporal groupings. After each new parturition, once we know the exact birthing time, the classifiers are re-evaluated to see which ones score best. Finally, the new data will be added as training data for the on-line classifier. 

## Video Analysis
Annotated video files, perhaps notebook with analysis to support dissociation-hypothesis. 






## Fiber-Photometry-Analysis




### Recording Analysis 




### Fiber-Position-Evaluation




## Brain-Behavior-Correlate-Analysis

