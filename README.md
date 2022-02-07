# Binaural auditory model RAA
Binaural auditory model using the AMT Toolbox framework, as used in the paper "Predicting perceived reverberation in diferent room acoustic environments using a binaural auditory model".

After you have installed this toolbox (see below) you just need to run the script **exp_osses2017.m**. This script will load the waveforms that are stored under the **auxdata** which in turn will be feed into the binaural model. The script **Run_RAA** will take the sound file names and will reshape/preprocess the sounds before they are fed into the binaural model, which is coded in the script **dorp2013.m** following the conventions of AMT 1.0.

# Installation
The following are the general instructions to get the binaural auditory model (RAA) for MATLAB operative in your computer. The toolbox has been tested on Windows and Linux, using MATLAB (versions R2012b-R2020b).

1. Download or clone the binaural-auditory-model-RAA project to your local computer (one way: press the button 'Code'->Choose 'Download ZIP' and unzip somewhere).
2. This toolbox requires the Auditory Modelling Toolbox v.1.0 (AMT 1.0) that can be downloaded from [here](http://amtoolbox.org/download.php). After the download you are not expected to do anything else, as the AMT toolbox will automatically be initialised in our next step:
3. Open and run the script **raastart.m**. This script will add all the paths under the fastACI toolbox to your local MATLAB path and it will run the script **amt_start.m** to initilise the AMT toolbox. 

# Citing this repository
This repository can be cited as follows:

Osses Vecchi, Alejandro (2017). "Binaural auditory model RAA," [https://github.com/aosses-tue/binaural-auditory-model-RAA/](https://github.com/aosses-tue/binaural-auditory-model-RAA/).
[![DOI](https://zenodo.org/badge/79450845.svg)](https://zenodo.org/badge/latestdoi/79450845)

# Other references
P. Majdak, C. Hollomey, & R. Baumgartner (2021). **AMT 1.0: The toolbox for reproducible research in auditory modeling**, submitted to Acta Acustica.

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
