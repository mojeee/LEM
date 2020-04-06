#include "flameSolver.h"
#include "scalarFunction.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>

using namespace std;

FlameSolver::FlameSolver()
    : U(0, 0, Stride1X(1 ,1))
    , T(0, 0, Stride1X(1 ,1))
    , Y(0, 0, 0, StrideXX(1 ,1))
    , jCorrSolver(jCorrSystem)
    , strainfunc(NULL)
    , rateMultiplierFunction(NULL)
    , stateWriter(NULL)
    , timeseriesWriter(NULL)
    , heatLossFunction(NULL)
    , tbbTaskSched(tbb::task_scheduler_init::deferred)
    , vzInterp(new BilinearInterpolator)
    , vrInterp(new BilinearInterpolator)
    , TInterp(new BilinearInterpolator)
{
}

FlameSolver::~FlameSolver()
{
    delete strainfunc;
    delete rateMultiplierFunction;
}

void FlameSolver::setOptions(const ConfigOptions& _options)
{
    SplitSolver::setOptions(_options);
    tStart = options.tStart;
    tEnd = options.tEnd;

    gas.setOptions(_options);
    grid.setOptions(_options);
}

void FlameSolver::initialize(void)
{
    tbbTaskSched.initialize (options.nThreads);
    delete strainfunc;
    strainfunc = newScalarFunction(options.strainFunctionType, options);

    delete rateMultiplierFunction;
    if (options.rateMultiplierFunctionType != "") {
        rateMultiplierFunction = newScalarFunction(options.rateMultiplierFunctionType,
                                                   options);
    } else {
        rateMultiplierFunction = NULL;
    }

    flamePosIntegralError = 0;
    terminationCondition = 1e10;

    // Cantera initialization
    gas.initialize();
    nSpec = gas.nSpec;
    nVars = nSpec + 2;
    W.resize(nSpec);
    gas.getMolecularWeights(W);

    // Get Initial Conditions
    loadProfile();
	
    grid.setSize(x.size());
    convectionSystem.setGas(gas);
    convectionSystem.setLeftBC(Tleft, Yleft);
    convectionSystem.setTolerances(options);

    for (size_t k=0; k<nVars; k++) {
        DiffusionSystem* term = new DiffusionSystem();
        TridiagonalIntegrator* integrator = new TridiagonalIntegrator(*term);
        integrator->resize(nPoints);
        diffusionTerms.push_back(term);
        diffusionSolvers.push_back(integrator);
    }
    if (options.wallFlux) {
        diffusionTerms[kEnergy].yInf = options.Tinf;
        diffusionTerms[kEnergy].wallConst = options.Kwall;
    }

    resizeAuxiliary();

    ddtConv.setZero();
    ddtDiff.setZero();
    ddtProd.setZero();
logFile.write(format("Hi moj 4444444444444444444444444444444444444444444444444"));
    updateChemicalProperties();
    calculateQdot();

    t = tStart;
    tOutput = t;
    tRegrid = t + options.regridTimeInterval;
    tProfile = t + options.profileTimeInterval;
    nTotal = 0;
    nRegrid = 0;
    nOutput = 0;
    nProfile = 0;
    nTerminate = 0;
    nCurrentState = 0;

    grid.updateValues();
    resizeAuxiliary();

    tFlamePrev = t;
    tNow = t;

    totalTimer.start();
}

void FlameSolver::setupStep()
{
    setupTimer.start();

    // Debug sanity check
    #ifndef NDEBUG
        bool error = false;
        for (size_t j=0; j<nPoints; j++) {
            if (T(j) < 295 || T(j) > 3000) {
                logFile.write(format(
                    "WARNING: Unexpected Temperature: T = %f at j = %i") % T(j) % j);
                error = true;
            }
        }
        if (error) {
            writeStateFile("err_setupStep", true, false);
        }
    #endif

    // Reset boundary conditions to prevent numerical drift
    if (grid.leftBC == BoundaryCondition::FixedValue) {
        T(0) = Tleft;
        Y.col(0) = Yleft;
    }

    if (grid.rightBC == BoundaryCondition::FixedValue) {
        T(jj) = Tright;
        Y.col(jj) = Yright;
    }
    logFile.write(format("Hi moj 3333333333333333333333333333333333333333333333333333"));
    updateChemicalProperties();

    updateBC();
    if (options.xFlameControl) {
        update_xStag(t, true); // calculate the value of rVzero
    }
    //convectionSystem.set_rVzero(rVzero);
    setupTimer.stop();

    // Set up solvers for split integration
    updateCrossTerms();
}

void FlameSolver::prepareIntegrators()
{
    splitTimer.resume();
    // Diffusion terms
    if (!options.quasi2d) {
        // Diffusion solvers: Energy and momentum
        diffusionTerms[kMomentum].B = rho.inverse();
        diffusionTerms[kEnergy].B = (rho * cp).inverse();

        diffusionTerms[kMomentum].D = mu;
        diffusionTerms[kEnergy].D = lambda;

        // Diffusion solvers: Species
        for (size_t k=0; k<nSpec; k++) {
            diffusionTerms[kSpecies+k].B =rho.inverse();
            diffusionTerms[kSpecies+k].D = rhoD.row(k);
        }
    } else {
        // Diffusion solvers: Energy and momentum
        diffusionTerms[kMomentum].B.setZero(nPoints);
        diffusionTerms[kEnergy].B.setZero(nPoints);

        diffusionTerms[kMomentum].D.setZero(nPoints);
        diffusionTerms[kEnergy].D.setZero(nPoints);

        // Diffusion solvers: Species
        for (size_t k=0; k<nSpec; k++) {
            DiffusionSystem& sys = diffusionTerms[kSpecies+k];
            sys.D = rhoD.row(k);
            for (size_t j = 0; j <= jj; j++) {
                sys.B[j] = 1 / (rho[j] * vzInterp->get(x[j], tNow));
            }
        }
    }

    setDiffusionSolverState(tNow);
    for (size_t i=0; i<nVars; i++) {
        diffusionTerms[i].splitConst = splitConstDiff.row(i);
    }

    // Production terms
    setProductionSolverState(tNow);
    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].splitConst = splitConstProd.col(j);
    }

    // Convection terms
    setConvectionSolverState(tNow);
    dmatrix ddt = ddtConv*0.0 + ddtDiff*0.0 + ddtProd;
    if (options.splittingMethod == "balanced") {
        ddt += ddtCross;
    }

    dvec tmp = (W.inverse().matrix().transpose() * ddt.bottomRows(nSpec).matrix()).array();
    drhodt = - rho * (ddt.row(kEnergy).transpose() / T + tmp * Wmx);

    assert(mathUtils::notnan(drhodt));
    //convectionSystem.setDensityDerivative(drhodt);
    //convectionSystem.setSplitConstants(splitConstConv);
    //convectionSystem.updateContinuityBoundaryCondition(qDot, options.continuityBC);
    splitTimer.stop();
}

int FlameSolver::finishStep()
{
    logFile.verboseWrite("done!");

    // *** End of Strang-split integration step ***
    correctMassFractions();

    t = tNow + dt;
    tNow += dt;

    nOutput++;
    nRegrid++;
    nProfile++;
    nTerminate++;
    nCurrentState++;

    if (debugParameters::debugTimesteps) {
        int nSteps = convectionSystem.getNumSteps();
        logFile.write(format("t = %8.6f (dt = %9.3e) [C: %i]") % t % dt % nSteps);
    }

    setupTimer.resume();
    if (t + 0.5 * dt > tOutput || nOutput >= options.outputStepInterval) {
        calculateQdot();
        timeVector.push_back(t);
        heatReleaseRate.push_back(getHeatReleaseRate());
        saveTimeSeriesData("out", false);
        tOutput = t + options.outputTimeInterval;
        nOutput = 0;
    }

    // Periodic check for terminating the integration (based on steady heat
    // release rate, etc.) Quit now to skip grid adaptation on the last step
    if (t >= tEnd) {
        return 1;
    } else if (nTerminate >= options.terminateStepInterval) {
        nTerminate = 0;
        if (checkTerminationCondition()) {
            return 1;
        }
    }




//**--------------------------------------- Main Part Of Code -----------------------------------------------

	if (init_flag==1)
	{
		INIT_AllParameters();

		init_flag=0;
	}
			
		double* F = new double[nPoints * nSpec];
		double* FE = new double[nPoints];
		double RHOM,RHOP,SUMYK,SUMX,TDOT,XMDXM;
		int k,j,kk;
	logFile.write(format("Hi moj 111111111111111111111111111111111111111111111111111111"));
	updateChemicalProperties();
	double** YVelocity;
	YVelocity=DiffusionVelocityCalculator();
	logFile.write(format("Hi moj, the counter number in the loop for solving flame is : %i ") %Check_flag );
			// PREMIXADV should be located here ------- PREMIX ADVANCEMENT IN TIME

// PREMIXADV FUNCTION -----------------------------------------------------------------------------------


//----------------------INTERIOR MESH POINTS-----------------------------

      //INTERIOR CELLS
      for( j = 1;j<nPoints;j++)
	{
               	
	        for (k = 0; k<nSpec; k++)
		{
	        	wDot(k,j) = wDot(k,j)*GFAC;
		}
		//-------------------------------SPECIES CONSERVATION EQUATION--------------------------------

	        SUMYK = 0.0;
		
		// set rho
		RHOP= rho(j);
		RHOM=rho(j-1);

		// XMDOT is set in SetIC
	        XMDXM = XMDOT / DX;
	        	for (k = 0; k<nSpec;k++)
	        	{
	           		SUMYK = SUMYK + Y(j,k);
				//species molecular weights === XMWT
	           		*(F + j*nSpec + k) = XMDXM * ( Y(k,j)-Y(k,j-1) )*0.0;// convection term set to zero
				*(F + j*nSpec + k) = *(F + j*nSpec + k) + (RHOP*YVelocity[j][k] - RHOM*YVelocity[j-1][k])/DX;
				*(F + j*nSpec + k) = *(F + j*nSpec + k) - wDot(k,j)*W(k);
				*(F + j*nSpec + k) = - dt*( *(F + j*nSpec + k) )/((RHOP + RHOM)/2.0);
	        	}



	}


	for (j = 1;j<nPoints;j++)
	 {
      		for(k = 0;k<nSpec;k++)
	 	 {
      			Y(k,j) = Y(k,j) + *(F + j*nSpec + k);
			if(Y(k,j)< 0.0 ) 
			 {
				Y(k,j) = 0.0;
			 }
	         }
         }


	delete[] F;

                //-------------------------------ENERGY EQUATION-----------------------------------------------
	for( j = 1;j<nPoints;j++)
	{
	       	for (k = 0; k<nSpec; k++)
		{
	        	wDot(k,j) = wDot(k,j)*GFAC;
		}

		SUMX = 0.0;
	        TDOT = 0.0;

		// set rho
		RHOP= rho(j);
		RHOM=rho(j-1);

	        for (k = 0; k<nSpec;k++)
		{
			TDOT = TDOT + wDot(k,j)*hk(k,j);
			SUMX = SUMX + 0.25 * (RHOP*YVelocity[j][k] + RHOM*YVelocity[j-1][k]) *cpSpec(k,j)*(T(j+1)-T(j-1))/DX;
		}
		


	        FE[j] = (XMDOT*(T(j)-T(j-1))/DX)*0.0;// convection term set to zero
		if ( j == ( nPoints - 1 ) )
		 {
			FE[j] = FE[j] +(lambda(j-1)*(T(j)-T(j-1))/DX)/(cp(j)*DX);
		 }
		else
		 {
			FE[j] = FE[j] -(lambda(j)*(T(j+1)-T(j))/DX-lambda(j-1)*(T(j)-T(j-1))/DX)/(cp(j)*DX);
		 }		
		FE[j] = FE[j] +(SUMX+TDOT)/cp(j);
		FE[j] = - dt*FE[j]/((RHOP+ RHOM)/2.0);

	}

	
	//----------------------- UPDATE ARRAYS --------------------------

	for (j = 1;j<nPoints;j++)
	 {
      		T(j) = T(j) + FE[j];

         }

	delete[] FE;
	delete[] YVelocity;
	
	 logFile.write(format("Hi moj, the NTS_PE is : %i ") %NTS_PE);

	if (Check_flag%200==0)
		{
			TM();
		}

	if (Check_flag%NTSPSIM==0 && NSIM_counter<=NSIM)
	{
		check_velocity=XMDOT*(8.3144627*T(nPoints-1))/(P*Wmx(nPoints-1));
        	if( check_velocity > 2.4 || check_velocity <= 0.0 )
		 {
			init_flag=1;
			error_flag=error_flag+1;	
        	 }
		else
		 {		
			Print_flag = 1;
			NSIM_counter=NSIM_counter+1;
		 }
	}
	
	Check_flag=Check_flag+1;







/*
// PREMIXADV FUNCTION -----------------------------------------------------------------------------------

			
			// triplet map			
			if (NTS_COUNT%NTS_PE == 0)
			{
				
				logFile.write(format("Hi moj, tripletmap is done ."));
			}
			
			//CHECK/ADJUST FUEL FLOW CFUEL Function			
			
			// tripletmap counter
			logFile.write(format("Hi moj, the triplet map counter is : %i ") %NTS_COUNT );
			NTS_COUNT=NTS_COUNT+1;
		}
		
	}
*/
/*	
	// PremixADV should be located here ------- PREMIX ADVANCEMENT IN TIME

	
	
	// CFUel function ---  CHECK/ADJUST FUEL FLOW

	

	// check for error in last cell and prepare for print
	if (Check_flag%NTSPSIM==0 && Print_counter<=NSIM)
	{
		check_velocity=XMDOT*(8.3144627*T(nPoints-1))/(P*Wmx(nPoints-1));
        	if( check_velocity > 2.4 || check_velocity <= 0.0 )
		 {
			init_flag=1;
			error_flag=error_flag+1;	
        	 }
		else
		 {		
			Print_flag = 1;
			Print_counter = Print_counter+1;
		 }
	}

	// show the error flag value in logfile
        logFile.write(format("Hi moj, the Error Flag value is :   %i") %error_flag );

	if ( Print_flag == 1 )
	{
		XRecord();
	}
*/
/*	if (t==2e-9)
	{
		INIT_AllParameters();
	}	

	double** YV;
	YV=DiffusionVelocityCalculator();
	
	if (t== 2e-9){
        ofstream debug_MassF ("debug_mass_fraction.csv");
	for (int j=0; j<nPoints;j++)
	{
        	
		for (int k=0; k<nSpec;k++)
		{
		
				debug_MassF << YV[j][k] << "\t";
		
		}
		
	}
		debug_MassF.close();
	}

    if (t > tRegrid || nRegrid >= options.regridStepInterval) {
        if (debugParameters::debugAdapt || debugParameters::debugRegrid) {
            writeStateFile("preAdapt", false, false);
        }
        regridTimer.start();
        tRegrid = t + options.regridTimeInterval;
        nRegrid = 0;

        // If the left grid point moves, a new boundary value for rVzero
        // needs to be calculated from the mass flux V on the current grid points
        dvec x_prev = grid.x;
        convectionSystem.evaluate();
        dvec V_prev = convectionSystem.V;

        // dampVal sets a limit on the maximum grid size
        grid.dampVal.resize(grid.x.rows());
        for (size_t j=0; j<nPoints; j++) {
            double num = std::min(mu[j],lambda[j]/cp[j]);
            for (size_t k=0; k<nSpec; k++) {
                if (rhoD(k,j) > 0) {
                    num = std::min(num, rhoD(k,j));
                }
            }
            double den = std::max(std::abs(rho[j]*strainfunc->a(t)), 1e-100);
            grid.dampVal[j] = sqrt(num/den);
        }
        dvec dampVal_prev = grid.dampVal;

        vector<dvector> currentSolution;
        rollVectorVector(currentSolution, state);
        rollVectorVector(currentSolution, ddtConv);
        rollVectorVector(currentSolution, ddtDiff);
        rollVectorVector(currentSolution, ddtProd);

        grid.nAdapt = nVars;
        if (options.quasi2d) {
            // do not change grid extents in this case
        } else if (strainfunc->a(tNow) == 0) {
            calculateQdot();
            grid.regridUnstrained(currentSolution, qDot);
        } else {
            grid.regrid(currentSolution);
        }

        // Interpolate dampVal onto the modified grid
        grid.dampVal = mathUtils::interp1(x_prev, dampVal_prev, grid.x, false);

        grid.adapt(currentSolution);

        // Perform updates that are necessary if the grid has changed
        if (grid.updated) {
            logFile.write(format("Grid size: %i points.") % nPoints);
            resizeMappedArrays();
            unrollVectorVector(currentSolution, state, 0);
            unrollVectorVector(currentSolution, ddtConv, 1);
            unrollVectorVector(currentSolution, ddtDiff, 2);
            unrollVectorVector(currentSolution, ddtProd, 3);
            correctMassFractions();

            // Update the mass flux (including the left boundary value)
            rVzero = mathUtils::interp1(x_prev, V_prev, grid.x[0]);
            convectionSystem.utwSystem.V = mathUtils::interp1(x_prev, V_prev, grid.x);

            // Allocate the solvers and arrays for auxiliary variables
            resizeAuxiliary();
            calculateQdot();

            if (debugParameters::debugAdapt || debugParameters::debugRegrid) {
                writeStateFile("postAdapt", false, false);
            }
            grid.updated = false;
        }
        regridTimer.stop();
    }
*/

    // Save the current integral and profile data in files that are
    // automatically overwritten, and save the time-series data (out.h5)
    if (nCurrentState >= options.currentStateStepInterval) {
        calculateQdot();
        nCurrentState = 0;
        saveTimeSeriesData("out", true);
        writeStateFile("profNow");
    }

    // *** Save flame profiles
    if (t + 0.5 * dt > tProfile || nProfile >= options.profileStepInterval) {
        if (options.outputProfiles) {
            writeStateFile();
        }
        tProfile = t + options.profileTimeInterval;
        nProfile = 0;
    }
    
setupTimer.stop();
    if (nTotal % 10 == 0) {
        printPerformanceStats();
    }
    nTotal++;
    return 0;
}


void FlameSolver::XRecord()

{
            int profile,j;
            profile=Print_counter;
	    string profileString=std::to_string(profile);
            string pref="result/prof00";
            string suf=".txt";
            string filename1= profileString+suf;
            string filename= pref+filename1;
            ofstream prof (filename);
	    prof<< "Temperature profile at time : " << "\t" << t << "\n" ;
       	    for ( j= 0; j< nPoints; j++)
        	 {
			prof<< T(j) << "\n" ;
         	}
            prof.close();
}

void FlameSolver::ReadParameters(Config& config)
 {
    ifstream fin("config.txt");
    string line;
    while (getline(fin, line)) {
        istringstream sin(line.substr(line.find("=") + 1));
        if (line.find("endtime") != -1)
            sin >> config.endtime;
        else if (line.find("timestep") != -1)
            sin >> config.timestep;
        else if (line.find("Re_t") != -1)
            sin >> config.Re_t;
        else if (line.find("dom") != -1)
            sin >> config.dom;
        else if (line.find("pressure") != -1)
            sin >> config.pressure;
        else if (line.find("u") != -1)
            sin >> config.velocity;
        else if (line.find("T") != -1)
            sin >> config.Temp; // GET_RHO_U
        else if (line.find("kinematic_viscosity") != -1)
            sin >> config.kinematic_viscosity;
        else if (line.find("GFAC") != -1)
            sin >> config.GFAC;
        else if (line.find("FAL") != -1)
	    sin >> config.FAL;
        else if (line.find("Intlength") != -1)
            sin >> config.Intlength;
        else if (line.find("NofRperR") != -1)
            sin >> config.NofRperR;
        else if (line.find("NSPE") != -1)
            sin >> config.NSPE;

    }
}
double** FlameSolver::DiffusionVelocityCalculator()
{
	int j,k;
	double SUM,VC;
        double Yp[nSpec][nPoints];

		
		for(j=0;j<nPoints-1;j++)
		{
			for(k=0;k<nSpec;k++)
			{			
				Yp[k][j] = Y(k,j+1);
			}
		}


		for(k=0;k<nSpec;k++)
		{			
			Yp[k][nPoints-1] = Y(k,nPoints-1);
		}


	double** YV= new double*[nPoints];

        for(j=0; j<nPoints; j++)
	{

		YV[j]= new double[nSpec];
		for(k=0;k<nSpec;k++)
		{
                        YV[j][k] = - Dkm(k,j)*(Yp[k][j]-Y(k,j))/DX;
		}

		SUM = 0.0;
		for(k=0;k<nSpec;k++)
		{
                        SUM = SUM + YV[j][k];
		}
                VC = - SUM;
		for(k=0;k<nSpec;k++)
		{
                        YV[j][k] = YV[j][k] + Y(k,j)*VC;
			
		}
	}
	
	if(t==2e-9){
        ofstream proof ("diff_velocity.txt");
	    
       	    for ( k= 0; k< nPoints; k++)
        	{
			proof<< YV[k][9] << "\n" ;
         	}
            proof.close();}
	return YV;
}

void FlameSolver::INIT_AllParameters() 
{
	Config config;
        ReadParameters(config);
	double Temperature,U_velocity,density;
	dt=config.timestep;
	NTS  = config.endtime/dt;	
	DOM = config.dom;
	NC= nPoints;
	NCM1=NC-1;
	NCP1=nPoints+1;
	DX=DOM/NCM1;
	XMDT = dt/DX;	
	GFAC=config.GFAC;
	NFL=config.FAL;
	NSIM =config.NofRperR ;
	NTSPSIM=config.NSPE ;
 	NTS_COUNT = 0;
	Temperature= config.Temp;
	U_velocity= config.velocity;
	P =config.pressure;
        density = P*Wmx(nPoints-1)/(8.3144627*Temperature);
        XMDOT    = density*U_velocity;
	XNU = config.kinematic_viscosity;
	Re =  config.Re_t;
	XLint = config.Intlength;
	XLk   = XLint/pow(Re,0.75);
	C_lambda = 15.0 ;
	Rate = DOM*100*(54.0/5.0)*( XNU*Re / (C_lambda*pow(XLint,3)) )*( pow((XLint/XLk),(5/3)) - 1)/( 1 - pow((XLk/XLint),(4/3)) );
	NTS_PE = (1.0/Rate)/dt+1;
	PDFA = pow(XLint,(5.0/3.0)) * pow(XLk,(-5.0/3.0)) / ( pow((XLint/XLk),(5.0/3.0)) -1.0 );
        PDFB = -pow(XLint,(5.0/3.0))/(  pow((XLint/XLk),(5.0/3.0)) -1.0 );
	TAU = XLk*XLk/XNU;
	if( (1.0/Rate)/dt < 1.0 )
	{
         MTS = int(Rate*dt)+1;
	}
	else
	{
         MTS = 1;
        }
	
}


void FlameSolver::Random_Number()
{
	 
	random=double(rand())/RAND_MAX;
	 
}



void FlameSolver::eddyLength()
{
        // generate a random number between 0 and 1
        Random_Number();

	int NSize;
        // make sure eddy is Greater than 6 cells long
        NSize = int(pow((random-PDFA)/PDFB,(-3.0/5.0))/DX);

        while (NSize<5)
        {
		Random_Number();
                NSize = int(pow((random-PDFA)/PDFB,(-3.0/5.0))/DX);
        	//logFile.write(format("hi moj, the size value   NSize=%i  random=%d") % NSize % random);

	}
        // make sure eddy size is divisible by 3
        if((NSize%3)==0)
        {
                L=NSize;
        }
        else if ((NSize%3)==1)
        {
                L=NSize-1;
        }
        else if ((NSize%3)==2)
        {
                L=NSize+1;
        }

}




void FlameSolver::BTriplet(double var[])
{
        // Permute Cells M through M+L-1 of the array S
        // as prescribrd by the discrete triplet map , where L is an integer multiple by 3.
        int Lo,k,j;
        Lo=int(L/3);
        double X[L];

        // first part of mapping
        for(j=1;j<=Lo;j++)
         {
                k=M+3*(j-1);
                X[j]=var[k]; //gather the cells going to the 1st image 
         }

        // second part of mapping
        for (j=1;j<=Lo;j++)
         {
                k=M+L+1-(3*j); // minus sign because second image is flipped
                X[j+Lo]=var[k]; // gather the cells going to 2nd image
         }

        // third part of mapping
        for (j=1;j<=Lo;j++)
         {
                k=M+(3*j)-1;
                X[j+Lo+Lo]=var[k];// gather the cells going to the 3rd image 
         }

        for(j=1;j<=L;j++)
         {
	        k=M+j-1;
	        var[k]=X[j];
         }

}





void FlameSolver::TM()
{
	// MTS shows number of triplet map require in each realization 
	// first of all call a random number 
	double Temp[nPoints];
        double y_x[nPoints];
        int j,k;
	
	// first of all call a random number
	Random_Number();
	
	// m determine the starting point of triplet map
	M=int(random*nPoints);

	// this loop check the starting point of triplet map
	while(M<(NCP1/4))
       	 {
		Random_Number();
		M=int(random*nPoints);
	 }
        
	// calculate the eddy length
	eddyLength();	

	// check eddy size does not exceed domain
	if((L+M)>nPoints)
         {
		M=nPoints-L;
	 }
	
	while((L+M)>nPoints || M<(NCP1/4))
	 {
		Random_Number();
		M=int(random*nPoints);
		eddyLength();
	 }
		
	for(j=0;j<nPoints;j++)
         {
                Temp[j]=T(j);
         }

         BTriplet(Temp);

        for(j=0;j<nPoints;j++)
         {
                T(j)=Temp[j];
         }

	for(k=0;k<=nSpec;k++)
         {
		for(j=0;j<nPoints;j++)
  		 {                        
			y_x[j]=Y(k,j);
		
		 }
		
		BTriplet(y_x);

		for(j=0;j<nPoints;j++)
		 {                        
			Y(k,j)=y_x[j];
		 }
	 }


	
}


void FlameSolver::PREMIXADV(double **YV)
{
	dmatrix F;	
	double RHOM,RHOP,SUMYK,SUMX,TDOT,XMDXM;
	int k,j,kk,jj;
        //EVALUATE AND STORE THE DIFFUSION VELOCITIES

			//double** YV;
			//YV=DiffusionVelocityCalculator();

      //----------------------INTERIOR MESH POINTS-----------------------------

      //INTERIOR CELLS
      for( j = 1;j<nPoints;j++)
	{
               	
	        for (k = 0; k<nSpec; k++)
		{
	        	wDot(k,j) = wDot(k,j)*GFAC;
		}

		//-------------------------------SPECIES CONSERVATION EQUATION--------------------------------

	        SUMYK = 0.0;
		
		// set rho
		RHOP= rho(j);
		RHOM=rho(j-1);

		// XMDOT is set in SetIC
	        XMDXM = XMDOT / DX;
	        	for (k = 0; k<nSpec;k++)
	        	{
	           		SUMYK = SUMYK + Y(j,k);
				//species molecular weights === XMWT
	           		F(3+k,j) = XMDXM * ( Y(k,j)-Y(k,j-1) )*0.0;// convection term set to zero
				F(3+k,j) = F(3+k,j) + (RHOP*YV[j][k] - RHOM*YV[j-1][k])/DX;
				F(3+k,j) = F(3+k,j) - wDot(k,j)*W(k);
				F(3+k,j) = - dt*F(3+k,j)/((RHOP + RHOM)/2.0);
	        	}

                //-------------------------------ENERGY EQUATION-----------------------------------------------

	        SUMX = 0.0;
	        TDOT = 0.0;

	        for (k = 0; k<nSpec;k++)
		{
			TDOT = TDOT + wDot(k,j)*hk(k,j);
			SUMX = SUMX + 0.25 * (RHOP*YV[j][k] + RHOM*YV[j-1][k]) *cpSpec(k,j)*(T(j+1)-T(j-1))/DX;
		}

	        F(3,j) = (XMDOT*(T(j)-T(j-1))/DX)*0.0;// convection term set to zero
		if ( j == ( nPoints - 1 ) )
		 {
			F(3,j) = F(3,j) +(lambda(j-1)*(T(j)-T(j-1))/DX)/(cp(j)*DX);
		 }
		else
		 {
			F(3,j) = F(3,j) -(lambda(j)*(T(j+1)-T(j))/DX-lambda(j-1)*(T(j)-T(j-1))/DX)/(cp(j)*DX);
		 }		
		F(3,j) = F(3,j) +(SUMX+TDOT)/cp(j);
		F(3,j) = - dt*F(3,j)/((RHOP+ RHOM)/2.0);

	}

	//----------------------- UPDATE ARRAYS --------------------------

	for (j = 1;j<nPoints;j++)
	 {
      		T(j) = T(j) + F(3,j);

      		for(k = 0;k<nSpec;k++)
	 	 {
      			Y(k,j) = Y(k,j) + F(3+k,j);
			if(Y(k,j)< 0.0 ) 
			 {
				Y(k,j) = 0.0;
			 }
	         }
         }

}

void FlameSolver::finalize()
{
    calculateQdot();
    saveTimeSeriesData("out", true);

    // *** Integration has reached the termination condition
    if (options.outputProfiles) {
        writeStateFile();
    }

    totalTimer.stop();
    printPerformanceStats();
    logFile.write(format("Runtime: %f seconds.") % totalTimer.getTime());
}

bool FlameSolver::checkTerminationCondition(void)
{
    if (tNow < options.tEndMin) {
        return false;
    }
    if (options.terminationMeasurement == "Q") {
        size_t j1 = mathUtils::findLast(timeVector < (tNow - options.terminationPeriod));
        if (j1 == npos) {
            logFile.write(format(
                    "Continuing integration: t (%8.6f) < terminationPeriod (%8.6f)") %
                    (tNow-timeVector[0]) % options.terminationPeriod);
            return false;
        }

        size_t j2 = timeVector.size()-1;
        double qMean = mathUtils::mean(heatReleaseRate,j1,j2);
        double hrrError = 0;
        for (size_t j=j1; j<=j2; j++) {
            hrrError += pow(heatReleaseRate[j]-qMean, 2);
        }
        hrrError = sqrt(hrrError) / (j2-j1+1);
        terminationCondition = hrrError / abs(qMean);

        logFile.write(format(
            "Heat release rate RMS error = %6.3f%%. absolute error: %9.4e") %
            (hrrError/qMean*100) % hrrError);

        if (hrrError/abs(qMean) < options.terminationTolerance) {
            logFile.write("Terminating integration: "
                    "Heat release RMS variation less than relative tolerance.");
            return true;
        } else if (hrrError < options.terminationAbsTol) {
            logFile.write("Terminating integration: "
                    "Heat release rate RMS variation less than absolute tolerance.");
            return true;
        }
    } else if (options.terminationMeasurement == "dTdt") {
        dvec dTdt = (ddtDiff + ddtConv + ddtProd + ddtCross).row(kEnergy);
        double value = (dTdt / T).matrix().norm() / sqrt(static_cast<double>(nPoints));
        logFile.write(format(
            "||1/T * dT/dt|| = %7.3f. Termination threshold = %7.2f") %
            value % options.termination_dTdtTol);
        if (value < options.termination_dTdtTol) {
            logFile.write("Terminating integration: "
                    "dTdt variation less than specified threshold.");
            return true;
        }
    }
    logFile.write(format(
        "Continuing integration. t = %8.6f") % (tNow-timeVector[0]));
    return false;
}

void FlameSolver::writeStateFile
(const std::string& fileNameStr, bool errorFile, bool updateDerivatives)
{
    if (stateWriter) {
        if (updateDerivatives) {
		logFile.write(format("Hi moj 55555555555555555555555555555555555555555555555555"));
            updateChemicalProperties();
            //convectionSystem.evaluate();
        }
        stateWriter->eval(fileNameStr, errorFile);
    }
}

void FlameSolver::saveTimeSeriesData(const std::string& filename, bool write)
{
    if (timeseriesWriter) {
        timeseriesWriter->eval(filename, write);
    }
}

void FlameSolver::resizeAuxiliary()
{
    size_t nPointsOld = rho.size();
    grid.setSize(T.size());

    if (nPoints == nPointsOld && !grid.updated) {
        return; // nothing to do
    }

    resizeMappedArrays();

    ddtCross.topRows(2).setZero();
    ddtCross.bottomRows(nSpec) *= NaN;

    rho.setZero(nPoints);
    drhodt.setZero(nPoints);
    Wmx.resize(nPoints);
    mu.resize(nPoints);
    lambda.setZero(nPoints);
    cp.setZero(nPoints);
    jCorr.resize(nPoints);
    sumcpj.setZero(nPoints);
    qDot.resize(nPoints);
    cpSpec.resize(nSpec, nPoints);
    rhoD.resize(nSpec, nPoints);
    Dkt.resize(nSpec, nPoints);
    Dkm.resize(nSpec, nPoints);
    wDot.resize(nSpec, nPoints);
    hk.resize(nSpec, nPoints);
    jFick.setZero(nSpec, nPoints);
    jSoret.setZero(nSpec, nPoints);

    grid.jj = nPoints-1;
    grid.updateBoundaryIndices();

    if (nPoints > nPointsOld) {
        for (size_t j=nPointsOld; j<nPoints; j++) {
            useCVODE.push_back(options.chemistryIntegrator == "cvode");

            SourceSystem* system;
            if (useCVODE[j]) {
                system = new SourceSystemCVODE();
            } else {
                system = new SourceSystemQSS();
            }
            // initialize the new SourceSystem
            system->setGas(&gas);
            system->initialize(nSpec);
            system->setOptions(options);
            system->setTimers(&reactionRatesTimer, &thermoTimer, &jacobianTimer);
            system->setRateMultiplierFunction(rateMultiplierFunction);
            system->setHeatLossFunction(heatLossFunction);
            system->setRhou(rhou);
            system->setPosition(j, x[j]);
            if (options.quasi2d) {
                system->setupQuasi2d(vzInterp, TInterp);
            }

            // Store the solver and system
            sourceTerms.push_back(system);
        }

    } else {
        // Delete the unwanted solvers and systems
        sourceTerms.erase(sourceTerms.begin()+nPoints, sourceTerms.end());
        useCVODE.erase(useCVODE.begin()+nPoints, useCVODE.end());
    }

    // Resize solution vector for diffusion systems / solvers
    for (size_t k=0; k<nVars; k++) {
        diffusionSolvers[k].resize(nPoints);
        diffusionTerms[k].setGrid(grid);
    }

    convectionSystem.setGrid(grid);
    convectionSystem.resize(nPoints, nSpec, state);
    convectionSystem.setLeftBC(Tleft, Yleft);

   // convectionSystem.utwSystem.setStrainFunction(strainfunc);
    //convectionSystem.utwSystem.setRhou(rhou);

    if (options.quasi2d) {
        convectionSystem.setupQuasi2D(vzInterp, vrInterp);
    }

    // Resize the jCorr stabilizer
    jCorrSolver.resize(nPoints);
    jCorrSystem.setGrid(grid);

    // Set the grid position for each of the source solvers
    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].setPosition(j, x[j]);
    }
}

void FlameSolver::resizeMappedArrays()
{
    resize(nVars, nPoints);
    remap(state, T, nPoints, kEnergy);
    remap(state, U, nPoints, kMomentum);
    remap(state, Y, nSpec, nPoints, kSpecies);
}

void FlameSolver::updateCrossTerms()
{
    assert(mathUtils::notnan(state));
    assert(mathUtils::notnan(rhoD));

    for (size_t j=0; j<jj; j++) {
        jCorr[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            jFick(k,j) = -0.5*(rhoD(k,j)+rhoD(k,j+1)) * ((Y(k,j+1)-Y(k,j))/hh[j]);
            jSoret(k,j) = -0.5*(Dkt(k,j)/T(j) + Dkt(k,j+1)/T(j+1))
                * (T(j+1)-T(j))/hh[j];
            jCorr[j] -= jFick(k,j) + jSoret(k,j);
        }
    }
    jCorr[jj] = 0;

    assert(mathUtils::notnan(jFick));
    assert(mathUtils::notnan(jSoret));
    assert(mathUtils::notnan(lambda));
    assert(mathUtils::notnan(rho));
    assert(mathUtils::notnan(cp));

    // Add a bit of artificial diffusion to jCorr to improve stability
    jCorrSolver.y = jCorr;
    jCorrSystem.B = lambda / (rho * cp); // treat as Le = 1
    jCorrSystem.D.setConstant(nPoints, 1.0);
    jCorrSystem.splitConst.setZero(nPoints);

    double dt = options.globalTimestep;
    for (size_t j=1; j<jj; j++) {
        dt = std::min(dt, options.diffusionTimestepMultiplier*dlj[j]*dlj[j]/(jCorrSystem.B[j]));
    }

    jCorrSolver.initialize(0, dt);
    jCorrSolver.integrateToTime(options.globalTimestep);
    assert(mathUtils::notnan(jCorrSolver.y));

    jCorr = jCorrSolver.y;

    // dYdt due to gradients in other species and temperature
    Eigen::Block<dmatrix> dYdtCross = ddtCross.bottomRows(nSpec);

    // dTdt due to gradients in species composition
    Eigen::Block<dmatrix, 1> dTdtCross = ddtCross.row(kEnergy);

    dYdtCross.col(0).setZero();
    dYdtCross.col(jj).setZero();
    for (size_t j=1; j<jj; j++) {
        sumcpj[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            dYdtCross(k,j) = -0.5 / (r[j] * rho[j] * dlj[j]) *
                (rphalf[j] * (Y(k,j) + Y(k,j+1)) * jCorr[j] -
                 rphalf[j-1] * (Y(k,j-1) + Y(k,j)) * jCorr[j-1]);
            dYdtCross(k,j) -= 1 / (r[j] * rho[j] * dlj[j]) *
                (rphalf[j] * jSoret(k,j) - rphalf[j-1] * jSoret(k,j-1));
            sumcpj[j] += 0.5*(cpSpec(k,j) + cpSpec(k,j+1)) / W[k] *
                (jFick(k,j) + jSoret(k,j) + 0.5 * (Y(k,j) + Y(k,j+1)) * jCorr[j]);
        }
        double dTdx = cfm[j] * T(j-1) + cf[j] * T(j) + cfp[j] * T(j+1);
        if (!options.quasi2d) {
            dTdtCross[j] = - 0.5 * (sumcpj[j] + sumcpj[j-1]) * dTdx / (cp[j] * rho[j]);
        }
    }

    assert(mathUtils::notnan(ddtCross));
    assert(mathUtils::notnan(sumcpj));
}

void FlameSolver::updateBC()
{
    BoundaryCondition::BC leftPrev = grid.leftBC;
    BoundaryCondition::BC rightPrev = grid.rightBC;

    if (options.wallFlux && x[0] >= 0.0 && x[0] <= options.centerGridMin) {
        grid.leftBC = BoundaryCondition::WallFlux;
    } else if (grid.ju == 0 || (grid.jb == 0 && grid.fixedBurnedVal)) {
        grid.leftBC = BoundaryCondition::FixedValue;
    } else if ((options.twinFlame || options.cylindricalFlame) &&
                x[0] >= 0.0 && x[0] <= options.centerGridMin) {
        grid.leftBC = BoundaryCondition::ControlVolume;
    } else {
        grid.leftBC = BoundaryCondition::ZeroGradient;
    }

    if (options.flameType == "premixed" && grid.jb == jj && !grid.fixedBurnedVal) {
        grid.rightBC = BoundaryCondition::Floating;
    } else {
        grid.rightBC = BoundaryCondition::FixedValue;
    }

    if (leftPrev != grid.leftBC) {
        logFile.write(format("updateBC: Left BC changed from %i to %i.") %
                      leftPrev % grid.leftBC);
    }

    if (rightPrev != grid.rightBC) {
        logFile.write(format("updateBC: Right BC changed from %i to %i.") %
                      rightPrev % grid.rightBC);
    }
}

void FlameSolver::updateChemicalProperties()
{
	logFile.write(format("Hi moj 66666666666666666666666666666666666666666666666"));
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nPoints,1),
                      TbbWrapper<FlameSolver>(&FlameSolver::updateChemicalProperties, this));
}

void FlameSolver::updateChemicalProperties(size_t j1, size_t j2)
{
    CanteraGas& gas = gases.local();
    if (!gas.initialized()) {
        gas.setOptions(options);
        gas.initialize();
    }
logFile.write(format("Hi moj 2"));
    // Calculate auxiliary data
    for (size_t j=j1; j<j2; j++) {
        // Thermodynamic properties
        thermoTimer.start();
        gas.setStateMass(&Y(0,j), T(j));
	
        rho[j] = gas.getDensity();
        Wmx[j] = gas.getMixtureMolecularWeight();
        cp[j] = gas.getSpecificHeatCapacity();
        gas.getSpecificHeatCapacities(&cpSpec(0,j));
        gas.getEnthalpies(&hk(0,j));
        thermoTimer.stop();

        // Transport Properties
        transportTimer.start();

        conductivityTimer.start();
        lambda[j] = gas.getThermalConductivity();
        conductivityTimer.stop();

        viscosityTimer.start();
        mu[j] = gas.getViscosity();
        viscosityTimer.stop();

        diffusivityTimer.start();
        gas.getWeightedDiffusionCoefficientsMass(&rhoD(0,j));
        gas.getThermalDiffusionCoefficients(&Dkt(0,j));
        gas.getDiffusionCoefficientsMole(&Dkm(0,j));
        diffusivityTimer.stop();
        transportTimer.stop();
    }
logFile.write(format("Hi moj 22222222222222222222222222222222222"));
}


void FlameSolver::setDiffusionSolverState(double tInitial)
{
    splitTimer.resume();
    size_t k = 0;
    for (TridiagonalIntegrator& integrator : diffusionSolvers) {
        double dt = (dlj.square().tail(nPoints - 2) /
            (diffusionTerms[k].B * diffusionTerms[k].D).segment(1, nPoints-2)).minCoeff();
        dt = std::min(options.diffusionTimestepMultiplier * dt, options.globalTimestep);
        integrator.initialize(tInitial, dt);
        k++;
    }

    for (size_t i=0; i<nVars; i++) {
        diffusionSolvers[i].y = state.row(i);
    }
    splitTimer.stop();
}

void FlameSolver::setConvectionSolverState(double tInitial)
{
    splitTimer.resume();
 //   convectionSystem.setState(tInitial);
    splitTimer.stop();
}

void FlameSolver::setProductionSolverState(double tInitial)
{
    splitTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].setState(tInitial, U(j), T(j), Y.col(j));
    }
    splitTimer.stop();
}

void FlameSolver::integrateConvectionTerms()
{
   /* setConvectionSolverState(tStageStart);
    convectionTimer.start();
    try {
        convectionSystem.integrateToTime(tStageEnd);
    } catch (DebugException& e) {
        logFile.write(e.errorString);
        writeStateFile("err_convectionIntegration", true, false);
        throw;
    }
    convectionTimer.stop();

    splitTimer.resume();
    convectionSystem.unroll_y();
    splitTimer.stop();*/
}

void FlameSolver::integrateProductionTerms()
{
    setProductionSolverState(tStageStart);
    reactionTimer.start();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nPoints,1),
                     TbbWrapper<FlameSolver>(&FlameSolver::integrateProductionTerms, this));
    reactionTimer.stop();
}

void FlameSolver::integrateProductionTerms(size_t j1, size_t j2)
{
    CanteraGas& gas = gases.local();
    if (!gas.initialized()) {
        gas.setOptions(options);
        gas.initialize();
    }

    int err = 0;
    for (size_t j=j1; j<j2; j++) {
        SourceSystem& system = sourceTerms[j];
        logFile.verboseWrite(format("%i") % j, false);
        system.setGas(&gas);
        if (int(j) == options.debugSourcePoint &&
            tStageEnd >= options.debugSourceTime) {
            system.setDebug(true);
            std::ofstream steps("sourceTermSteps.py");
            system.writeState(steps, true);

            while (system.time() < tStageEnd - tStageStart && err >= 0) {
                err = system.integrateOneStep(tStageEnd - tStageStart);
                system.writeState(steps, false);
            }

            system.writeJacobian(steps);
            steps.close();
            std::terminate();

        } else {
            err = system.integrateToTime(tStageEnd - tStageStart);
        }
        if (err >= 0) {
            logFile.verboseWrite(format(" [%s]...") % system.getStats(), false);
            sourceTerms[j].unroll_y();
            U(j) = sourceTerms[j].U;
            T(j) = sourceTerms[j].T;
            Y.col(j) = sourceTerms[j].Y;
        } else {
            // Print gas mole fractions to help identify problematic reactions
            dvec X0(nSpec), X1(nSpec);
            gas.thermo.setMassFractions_NoNorm(Y.col(j).data());
            gas.getMoleFractions(X0);
            gas.thermo.setMassFractions_NoNorm(sourceTerms[j].Y.data());
            gas.getMoleFractions(X1);
            logFile.write(format("Error at j = %i. Gas state:\n") % j);
            logFile.write(" k        Species     X Initial    X Final     Delta X");
            logFile.write("----  --------------  ----------  ----------  ----------");
            logFile.write(format("      %14s  %10.4g  %10.4g  %10.4g") %
                "T" % T(j) % sourceTerms[j].T % (sourceTerms[j].T - T(j)));
            for (size_t k = 0; k < nSpec; k++) {
                if (std::abs(X0[k]) > 1e-4 ||
                    std::abs(X1[k]) > 1e-4 ||
                    std::abs(X1[k]-X0[k]) > 1e-6) {
                    logFile.write(format("%4s  %14s  %10.4g  %10.4g  %10.4g") %
                        k % gas.thermo.speciesName(k) % X0[k] % X1[k] % (X1[k]-X0[k]));
                }
            }
            if (debugParameters::veryVerbose) {
                logFile.write(format("\nT = %s") % system.T);
                logFile.write(format("U = %s") % system.U);
                logFile.write("Y = ", false);
                Eigen::IOFormat fmt(15, Eigen::DontAlignCols, ", ", ", ", "", "", "[", "]");
                logFile.write(system.Y.format(fmt));
                if (options.nThreads == 1) {
                    // This diagnostic file can only be written when running with a
                    // single thread to avoid calling Python from threads that were
                    // not initialized to use Python
                    writeStateFile((format("prod%i_error_t%.6f_j%03i") %
                            1 % tStageEnd % j).str(), true, false);
                }
            }
        }
    }
}

void FlameSolver::integrateDiffusionTerms()
{
    setDiffusionSolverState(tStageStart);
    diffusionTimer.start();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nVars, 1),
                      TbbWrapper<FlameSolver>(&FlameSolver::integrateDiffusionTerms, this));
    diffusionTimer.stop();
}

void FlameSolver::integrateDiffusionTerms(size_t k1, size_t k2)
{
    for (size_t k=k1; k<k2; k++) {
        diffusionSolvers[k].integrateToTime(tStageEnd);
        assert(mathUtils::almostEqual(diffusionSolvers[k].t, tStageEnd));
        state.row(k) = diffusionSolvers[k].y;
    }
}

void FlameSolver::rollVectorVector
(vector<dvector>& vv, const dmatrix& M) const
{
    size_t N = vv.size();
    vv.resize(N + nVars, dvector(nPoints));

    for (size_t k=0; k<nVars; k++) {
        Eigen::Map<dvec>(&vv[N+k][0], nPoints) = M.row(k);
    }
}


void FlameSolver::unrollVectorVector
(vector<dvector>& vv, dmatrix& M, size_t i) const
{
    for (size_t k=0; k<nVars; k++) {
        M.row(k) = Eigen::Map<dvec>(&vv[i*nVars+k][0], nPoints);
    }
}


void FlameSolver::update_xStag(const double t, const bool updateIntError)
{
    calculateQdot();
    xFlameActual = getFlamePosition();
    xFlameTarget = targetFlamePosition(t);
    if (updateIntError) {
        flamePosIntegralError += (xFlameTarget-xFlameActual)*(t-tFlamePrev);
        tFlamePrev = t;
    }

    // controlSignal is approximately a*xStag
    double controlSignal = options.xFlameProportionalGain *
        ((xFlameTarget - xFlameActual) +
        (flamePosIntegralError + (xFlameTarget - xFlameActual) * (t - tFlamePrev)) *
        options.xFlameIntegralGain);

    if (debugParameters::debugFlameRadiusControl) {
        logFile.write(format(
            "rFlameControl: rF=%g;  control=%g;  P=%g;  I=%g;  dt=%g") %
            xFlameActual %
            controlSignal %
            (options.xFlameProportionalGain * (xFlameTarget - xFlameActual)) %
            (options.xFlameProportionalGain * flamePosIntegralError *
            options.xFlameIntegralGain) %
            (t - tFlamePrev));
    }

    double a = strainfunc->a(t);
    if (alpha == 1) {
        rVzero = 0.5*rhoLeft*(controlSignal*abs(controlSignal)-a*x[0]*x[0]);
    } else {
        rVzero = rhoLeft*(controlSignal-a*x[0]);
    }
}


double FlameSolver::targetFlamePosition(double t)
{
    if (t <= options.xFlameT0) {
        return options.xFlameInitial;
    } else if (t >= options.xFlameT0 + options.xFlameDt) {
        return options.xFlameFinal;
    } else {
        return options.xFlameInitial + (options.xFlameFinal - options.xFlameInitial) *
            (t - options.xFlameT0) / options.xFlameDt;
    }
}


void FlameSolver::calculateQdot()
{
    reactionRatesTimer.start();
    if (rateMultiplierFunction) {
        gas.setRateMultiplier(rateMultiplierFunction->a(tNow));
    }
    for (size_t j=0; j<nPoints; j++) {
        gas.setStateMass(&Y(0,j), T(j));
        gas.getEnthalpies(&hk(0,j));
        gas.getReactionRates(&wDot(0,j));
        qDot[j] = - (wDot.col(j) * hk.col(j)).sum();
    }
    reactionRatesTimer.stop();
}


void FlameSolver::correctMassFractions() {
    setupTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        gas.setStateMass(&Y(0,j), T(j));
        gas.getMassFractions(&Y(0,j));
    }
    setupTimer.stop();
}

double FlameSolver::getHeatReleaseRate(void)
{
    return mathUtils::integrate(x, qDot);
}

double FlameSolver::getConsumptionSpeed(void)
{
    double QoverCp = mathUtils::integrate(x, qDot / cp);
    double rhouDeltaT = rhou*(T(grid.jb)-T(grid.ju));
    return QoverCp/rhouDeltaT;
}

double FlameSolver::getFlamePosition(void)
{
    return mathUtils::trapz(x, (x*qDot).eval())/mathUtils::trapz(x,qDot);
}

void FlameSolver::loadProfile(void)
{
    grid.unburnedLeft = options.unburnedLeft;

    // Read initial condition specified in the configuration file
    x = options.x_initial;
    nPoints = x.size();
    resizeMappedArrays();

    U = options.U_initial;
    T = options.T_initial;
    if (options.Y_initial.rows() == static_cast<dmatrix::Index>(x.size())) {
        Y = options.Y_initial.transpose();
    } else {
        Y = options.Y_initial;
    }
    //convectionSystem.V = options.V_initial;
   // convectionSystem.utwSystem.V = options.V_initial;
    //rVzero = convectionSystem.utwSystem.V[0];

    grid.setSize(x.size());
    grid.updateValues();
    grid.updateBoundaryIndices();

    tNow = (options.haveTStart) ? options.tStart : 0.0;

    if (options.flameType == "premixed") {
        gas.setStateMass(&Y(0,grid.ju), T(grid.ju));
        rhou = gas.getDensity();
        gas.setStateMass(&Y(0,grid.jb), T(grid.ju));
        rhob = gas.getDensity();

        if (options.unburnedLeft) {
            rhoLeft = rhou;
            Tleft = T(grid.ju);
            Yleft = Y.col(grid.ju);
            rhoRight = rhob;
            Tright = T(grid.jb);
            Yright = Y.col(grid.jb);
        } else {
            rhoLeft = rhob;
            Tleft = T(grid.jb);
            Yleft = Y.col(grid.jb);
            rhoRight = rhou;
            Tright = T(grid.ju);
            Yright = Y.col(grid.ju);
        }

    } else if (options.flameType == "diffusion") {
        // Fuel composition
        size_t jFuel = (options.fuelLeft) ? 0 : jj;
        gas.thermo.setState_TPY(T(jFuel), options.pressure, &Y(0,jFuel));
        double rhoFuel = gas.getDensity();
        dvec Yfuel(nSpec);
        gas.getMassFractions(Yfuel);
        double Tfuel = gas.thermo.temperature();

        // Oxidizer composition
        size_t jOxidizer = (options.fuelLeft) ? jj : 0;
        gas.thermo.setState_TPY(T(jOxidizer),
                                options.pressure,
                                &Y(0,jOxidizer));
        double rhoOxidizer = gas.getDensity();
        dvec Yoxidizer(nSpec);
        gas.getMassFractions(Yoxidizer);
        double Toxidizer = gas.thermo.temperature();

        if (options.fuelLeft) {
            rhoLeft = rhoFuel;
            Tleft = Tfuel;
            Yleft = Yfuel;
            rhoRight = rhoOxidizer;
            Tright = Toxidizer;
            Yright = Yoxidizer;
        } else {
            rhoLeft = rhoOxidizer;
            Tleft = Toxidizer;
            Yleft = Yoxidizer;
            rhoRight = rhoFuel;
            Tright = Tfuel;
            Yright = Yfuel;
        }

        rhou = rhoOxidizer;

    } else if (options.flameType == "quasi2d") {
        gas.thermo.setState_TPY(T(0), options.pressure, &Y(0,0));
        rhoLeft = gas.thermo.density();
        Tleft = T(0);
        Yleft.resize(nSpec);
        gas.getMassFractions(Yleft);
        gas.thermo.setState_TPY(T(jj), options.pressure, &Y(0,jj));
        rhoRight = gas.thermo.density();
        Tright = T(jj);
        Yright.resize(nSpec);
        gas.getMassFractions(Yright);

        rhou = rhoRight;

    } else {
        throw DebugException("Invalid flameType: " + options.flameType);
    }

    updateBC();

    if (grid.leftBC == BoundaryCondition::ControlVolume &&
        options.xFlameControl)
    {
        double controlSignal;
        if (alpha == 0) {
            controlSignal = rVcenter/rhoLeft;
        } else {
            double tmp = pow(x[0],2) + 2*rVcenter/rhoLeft;
            controlSignal = mathUtils::sign(tmp)*sqrt(abs(tmp));
        }
        flamePosIntegralError = controlSignal /
            (options.xFlameProportionalGain * options.xFlameIntegralGain);
    }
}

void FlameSolver::printPerformanceStats(void)
{
    std::string filename = options.outputDir + "/stats";
    std::ofstream stats(filename.c_str(), std::ios::trunc | std::ios::out);
    totalTimer.stop();
    totalTimer.resume();
    stats << "\n   *** Performance Stats ***       time   ( call count )\n";
    printPerfString(stats, "                General Setup: ", setupTimer);
    printPerfString(stats, "             Split Term Setup: ", splitTimer);
    printPerfString(stats, "              Grid Adaptation: ", regridTimer);
    printPerfString(stats, "    Reaction Term Integration: ", reactionTimer);
    printPerfString(stats, "   Diffusion Term Integration: ", diffusionTimer);
    printPerfString(stats, "  Convection Term Integration: ", convectionTimer);
    printPerfString(stats, "                        Total: ", totalTimer);
    stats << "\n Subcomponents:\n";
    printPerfString(stats, "               Reaction Rates: ", reactionRatesTimer);
    printPerfString(stats, "         Transport Properties: ", transportTimer);
    printPerfString(stats, "          - thermal cond.    : ", conductivityTimer);
    printPerfString(stats, "          - viscosity        : ", viscosityTimer);
    printPerfString(stats, "          - diffusion coeff. : ", diffusivityTimer);
    printPerfString(stats, "     Thermodynamic Properties: ", thermoTimer);
    printPerfString(stats, "   Source Jacobian Evaluation: ", jacobianTimer);
    printPerfString(stats, "   UTW Convection Integration: ", convectionSystem.utwTimer);
    printPerfString(stats, "    Yk Convection Integration: ", convectionSystem.speciesTimer);
}

void FlameSolver::printPerfString(std::ostream& stats, const std::string& label,
                                  const PerfTimer& T)
{
    if (T.getTime() > 0.0) {
        stats << format("%s %9.3f (%12i)\n") % label % T.getTime() % T.getCallCount();
    }
}
