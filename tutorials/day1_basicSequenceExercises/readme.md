# Basic Pulseq Tutorial
## 1. Prerequisites
### 1.1 Programming Tools
Matlab software (https://mathworks.com/products/matlab.html) needs to be installed in your computer.   
Basic familiarity with Matlab programming is required.   
### 1.2 Pulseq Software
Please download the Pulseq software from https://github.com/pulseq/pulseq. Install Pulseq software in Matlab by adding its directory and subdirectories to Matlab's path.   
### 1.3 mapVBVD Software
The mapVBVD software is required to read raw data in Siemens TWIX format. The software can be downloaded from https://github.com/pehses/mapVBVD and installed by adding its directory to Matlab's path.
### 1.4 Text Compare Tool
We recommend using a text comparison tool to compare the sequences within the subsequent steps to visualize the changes that occur at each
step. *Meld* software can be used for text comparison and can be downloaded from <https://meldmerge.org/>. By using Meld or other comparable software packages, it is possible to very quickly make the relevant changes easily visible, as in the example below.

![meld](https://github.com/pulseq/ISMRM-Virtual-Meeting--November-15-17-2023/assets/26165904/306150db-68d7-4a8b-8eb3-13b8fccfc3a2)


## 2. Sequence Folder
### 2.1 Basic MR Spectroscopy
* *ex01_fid*: the simplest free induction decay (FID) sequence;   
* *ex02_fid2se*: spin-echo (SE) sequence without gradients;   
* *ex03_se_crushers*: SE sequence with a pair of crushers to eliminate spurious signals arising from the imperfect 180-degree RF pulse.
* *ex04_fid2gre1d*: from FID to 1D gradient recalled echo (GRE) sequence.   
### 2.2 Basic MR Imaging
* *ex11_gre1d2gre2d*: from 1D GRE to basic 2D GRE sequence;   
* *ex12_gre2d_lbl*: To implement labels to basic 2D GRE (ex11);   
* *ex13_gre2d_gradSpoil*: to implement gradient spoiling in readout and slice-selective directions for ex12;
* *ex14_gre2d_PErefocus*: to implement refocusing in phase-encoding direction for ex13;
* *ex15_gre2d_RFspoil*: to implement RF spoiling for ex14;
* *ex16_flash_fasterTming*: to shorten the duration of ex15 by doing "gradient surgery";
* *ex17_flash_accelerateComputation*: to accelerate the computation time for ex16 (optional). 

Note: *solXX* are the corresponding solution to the exercise *exXX*.   
## 3. Data and recon Folder
It contains the raw data, DICOM images, and reconstruction of each exercise.     

## 4. A More Detailed Pulseq Tutorial
This Pulseq tutorial only covers very basic sequence design concepts. For more detailed tutorials, please go to https://github.com/pulseq/tutorials.    
