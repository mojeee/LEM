#!/usr/bin/python
from ember import *

conf = Config(
    Paths(outputDir='Turbulence_K74'),
    InitialCondition(restartFile='prof000125.h5'),
    StrainParameters(initial=0,final=0),
    General(fixedLeftLocation = True,continuityBC='fixedLeft',flameGeometry='planar',nThreads=3),
    Chemistry(mechanismFile='Smooke.cti'),    
    Grid(vtol=0.1, dvtol=0.1,rmTol=0.0),
    TerminationCondition(measurement=None,tEnd = 0.05),
    Times(profileStepInterval=20,globalTimestep =2e-07,currentStateStepInterval = 10,regridTimeInterval=0.5,regridStepInterval=1000000),
    OutputFiles(heatReleaseRate= False, timeDerivatives= False, extraVariables=False, auxiliaryVariables=False),
    Debug(adaptation = False,regridding = False))
if __name__ == '__main__':
   
    conf.run()

