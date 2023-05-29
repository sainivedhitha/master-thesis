# master-thesis
Estimation of Joint Stiffness via a Musculoskeletal Model Driven by Motor Neuron Twitch Properties

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
