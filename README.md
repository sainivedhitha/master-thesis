# Estimation of Joint Stiffness via a Musculoskeletal Model Driven by Motor Neuron Twitch Properties

This research was undertaken as part of the MSc program. 

## Abstract 

Joint stiffness estimation involves joint perturbation or torque-based experiments in conjunction with system identification techniques. The modeling techniques recently used to estimate joint stiffness include bipolar ElectroMyoGraphy (EMG) data to drive the musculoskeletal model. These models do not provide detailed information on individual Alpha Motor Neuron ($\alpha$-MN) properties which is essential to improve the personalization of the Neuro Musculo Skeletal (NMS) model. The use of bipolar electrodes limits the resolution to extract the $\alpha$-MN properties. In this study, the NMS models were driven by the activation dynamics of EMG envelopes and motor unit activation dynamics, respectively. The normalized EMG envelopes were regarded as the activation profiles of EMG envelopes. For the activation dynamics of motor units, the Motor Units (MUs) were decomposed from uni-polar HD-EMG data using a blind source separation technique, and the motor unit distributions were sampled to formulate the activation profiles. The experimental torque from the Achilles and stiffness estimation from System Identification (SI) technique were used as a reference to validate the results of the models at the torque and stiffness level. The torque estimations by both models were better than the stiffness estimations. At the torque level, the model driven by motor units produced improved results; at the stiffness level, the model driven by EMG envelopes produced better results. Overall, this method can be used to predict the torque but enhancements should be made to increase the stiffness estimation results. Further, the entire study was performed under isometric conditions. The inclusion of unique properties of isometric conditions such as Short Range Stiffness (SRS) in the future could improve the results.

## Contents in this Repository

- PDF file containing the final version of the MSc thesis (scientific paper format).

- 'Result - Images' folder. This folder contains the figures of the results
	- Activation Dynamics (AD) of EMG Envelopes of all the 5 muscles per subject
	- AD of motor units per subject
	- Stiffness estimated by System Identification per subject
	- Torque estimated by the models per subject
	- Stiffness estimated by the models per subject
	- Average of all the subjects
	- Images of PNR, motor units, validation
	- pdf file containing the summary of the entire study

- 'Scripts for AD' folder. This folder contains the scripts required for calculating the AD of 
  EMG Envelopes and motor units.
	- 'parallelDecompositionCKC_MUinhibition.m' : Code that determines the AD of EMG Envelopes 
	  (normalized EMG Envelopes) and decomposes the motor units and finds parameters like discharge rate, recruitment threshold.
	  This code also perform quality control on the motor units to be used further. This code uses functions 'syncEMGandAchilles.m', 
	  'getSpikeTrains.m' and 'first_elements.m' along with the raw data and MVC data available in 'Experimental data' folder. The other 
	  functions required to run this code is in 'DEMUSE' folder.
	- 'Activation_Dynamics.m' : Code that determines PCA and the AD of the motor units. This uses the results of
	  the previous code.

- 'System Identification' folder. This folder contains the scripts required for calculating the joint stiffness using system identification.
	- 'TorquetoPRBS.m': Code which transforms the torque perturbation into a PRBS position signal.
	- 'MS_static.m' : Code which estimates stiffness from the raw experimental data. This code includes functions 'alignment.m', 'oulierRemoval.m'
	  and the functions of folder 'DL_tools'. This also uses the PRBS data obtained using the previous code.

- 'Stiffness.m' : Code that segments, normalizes and saves all the necessary data for calibration and execution of CEINMS.
	- This code uses the results in the folders 'Data\subject0x\AD' and 'Data\subject0x\SI'. The parameters segmented and saved includes the AD 
	  of EMG Envelopes of Tibialis Anterior, Gastrocnemius Medialis, Gastrocnemius Lateralis and Peroneus Longus, AD of motor units of Tibialis 
	  Anterior, Stiffness estimated by System Identification, experimental torque and experimental angle.
	- For calibration, the data were segmented, normalized and stored as 2 seconds data of length 1000.
	- For execution, the data were saved as entire 40 seconds data.

- 'Validation.m' : Code that loads the torque and stiffness results of CEINMS and validates them against the experimental torque and stiffness by 
  System Identification. The validation was performed using R-squared, RMSE and nRMSE values as well as through visual inspection of plots.
