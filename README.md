Turbulent Premixed Combustion modeling using The linear eddy model
Version 1.0
==================================================================
Author: Mojtaba Amini

Email: Mojtaba.amini.1995@gmail.com

Supervisor: Professor Mohammad Mahdi Salehi

Abstract
-------------------------


Direct calculation of production or consumption rate of all species in the reaction is a challenging task. Therefore the statistical method such as RANS or LES is used for calculating species changing rate. Also, some methods like FLAMELET would be used to simulate the interactions between kinetic and turbulence. Therefore the LEM (Linear Eddy Model) was implemented in the CANTERA (Open Source Library) to simulate one-dimensional and semi one-dimensional ﬂame with diﬀerent kinetic. The viscosity and molecular diﬀusion are implemented similar to a laminar ﬂame, and turbulence eﬀect is implemented by using Triplet Map (TM) (Stochastic advection). The LEM model can simulate a wide variety of ﬂames with different fuels and physics. In the end, the model creates a database like DNS, which uses the statistical results of the database to determine the PDF, SDR, wasting rate of variable, and comparing with DNS and experimental results


For each record, it is provided:

Files
-------------------------
The dataset includes the following files:

Files  | Description
------------- | -------------
'README.txt'  | 
'All_Matlab_Code.m'  | Almost all the essential function for post-processing is provided in each section. Also, there is enough comment for each section.
'COnditionalAverageCompare.m'  | For conditional Average calculation
'LaminarFlameSpeed.py'  | It calculates the flame speed with CANTERA code for the first estimation.
'F1_flame_case_for_Triplet.py'  | It will be used for initial condition calculation.
'InstructionToRun.txt'  | A short instruction to run LEM code
'K6_LEM.py'  | It will be used to run a job. 
'ModifiedreadDNS.m' | For reading DNS data
'Smooke.cti'  | Reduced mechanism that is used in the calculation.
'config.txt' | Problem configuration that is used to run a job by 'K6_LEM.py'
'flameSolver.cpp'  | Ember function that will be replaced with original files ( LEM functions are implemented here )
'flameSolver.h'  | Ember function that will be replaced with original files ( LEM functions are implemented here )



Notes
--------------------------
* Data will be captured in specific time, find more information in 'K6_LEM.py' and EMBER website.
       [EMBER Repository](http://speth.github.io/ember-doc/sphinx/html/index.html "EMBER Repository")
* This function jus uses the chemistry input of the EMBER code and don't use any function in EMBER
* To reduce the time it is necessary to developed the code separately
* For a problem with 3850 nodes and 18 species on a normal computer systems, it will take almost 2 days to get output.

Mojtaba Amini. January 2021.
