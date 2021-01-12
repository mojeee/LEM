#!/usr/bin/python
from ember import *

conf = Config(
    Paths(outputDir='DNSLaminarK6'),
    InitialCondition(flameType = 'premixed',pressure = 101325,fuel='CH4:1.0',oxidizer='O2:1, N2:3.76',
    nPoints = 1032,Tu = 300,equivalenceRatio=0.6,xLeft=0.0,xRight=0.01),
    StrainParameters(initial=0,final=0),
    General(fixedLeftLocation = True,continuityBC='fixedLeft',flameGeometry='planar',nThreads=3),
    Chemistry(mechanismFile='Smooke.cti'),
    Grid(vtol=0.1, dvtol=0.1,rmTol=0.0),
    TerminationCondition(tEnd = 0.0015),
    Times(profileStepInterval=200,globalTimestep = 6e-08,currentStateStepInterval = 100,regridTimeInterval=0.5,regridStepInterval=1000000),
    OutputFiles(heatReleaseRate= False, timeDerivatives= False, extraVariables=False, auxiliaryVariables=False),
    Debug(adaptation = False,regridding = False))
    

if __name__ == '__main__':
   
    conf.run()

