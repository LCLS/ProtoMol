#include <protomol/ProtoMolApp.h>
#include <protomol/base/Report.h>
#include <protomol/integrator/STSIntegrator.h>
#include <protomol/integrator/normal/NormalModeDiagonalize.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <algorithm>
#include <cmath>
#include <fstream>

#include <protomol/parallel/Parallel.h>
#include "protomol/base/Timer.h"
#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/integrator/Integrator.h"
#include "protomol/integrator/MTSIntegrator.h"
#include "protomol/integrator/StandardIntegrator.h"
#include "protomol/integrator/hessian/BlockHessian.h"
#include "protomol/integrator/hessian/BlockHessianDiagonalize.h"
#include "protomol/integrator/normal/NormalModeUtilities.h"
#include "protomol/type/EigenvectorInfo.h"
#include "protomol/type/SimpleTypes.h"

using namespace std;
using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol
{
  //__________________________________________________ NormalModeDiagonalize

  const string NormalModeDiagonalize::keyword( "NormalModeDiagonalize" );

  NormalModeDiagonalize::NormalModeDiagonalize() :
    MTSIntegrator(), NormalModeUtilities(), firstDiag(true),
    fullDiag(0), removeRand(0), rediagCount(0), nextRediag(0),
    validMaxEigv(0), myNextNormalMode(0), myLastNormalMode(0),
    rediagHysteresis(0), hessianCounter(0), rediagCounter(0),
    rediagUpdateCounter(0), eigenValueThresh(0), blockCutoffDistance(0),
    blockVectorCols(0), residuesPerBlock(0), memory_Hessian(0),
    memory_eigenvector(0), checkpointUpdate(false), origCEigVal(0),
    origTimestep(0), autoParmeters(false), adaptiveTimestep(0),
    postDiagonalizeMinimize(0), minLim(0), maxMinSteps(0), 
    geometricfdof(false), numerichessians(false) {
  }

  NormalModeDiagonalize::
  NormalModeDiagonalize(int cycles, int redi, bool fDiag, bool rRand,
                        Real redhy, Real eTh, int bvc, int rpb, Real dTh, 
                        bool apar, bool adts, bool pdm, Real ml, int maxit,
                        bool geo, bool num,
                        ForceGroup *overloadedForces,
                        StandardIntegrator *nextIntegrator ) :
    MTSIntegrator( cycles, overloadedForces, nextIntegrator ),
    NormalModeUtilities( 1, 1, 91.0, 1234, 300.0 ), firstDiag(true),
    fullDiag( fDiag ), removeRand( rRand ), rediagCount( redi ), nextRediag(0),
    validMaxEigv(0), myNextNormalMode(0), myLastNormalMode(0),
    rediagHysteresis( redhy ), hessianCounter( 0 ), rediagCounter( 0 ),
    rediagUpdateCounter(0), eigenValueThresh( eTh ),
    blockCutoffDistance( dTh ), blockVectorCols( bvc ),
    residuesPerBlock( rpb ),  memory_Hessian(0), memory_eigenvector(0),
    checkpointUpdate( false ), origCEigVal(0), origTimestep(0),
    autoParmeters(apar), adaptiveTimestep( adts ), postDiagonalizeMinimize(pdm),
    minLim(ml), maxMinSteps(maxit), geometricfdof(geo), numerichessians(num) {

    //find forces and parameters
    rHsn.findForces( overloadedForces );

  }

  NormalModeDiagonalize::~NormalModeDiagonalize()
  {
    //output stats
    report.precision( 5 );
    if ( rediagCounter && hessianCounter ) {
      report << plain
      << "NML Timing: Hessian: " << ( blockDiag.hessianTime.getTime() ).getRealTime() << "[s] (" << hessianCounter << " times)"
      << " diagonalize: " << ( blockDiag.rediagTime.getTime() ).getRealTime() << "[s] (" << rediagCounter << " re-diagonalizations)."
      << endl;

      if ( !fullDiag ) {
        report << plain << "NML Memory: "
        << "Hessian: "     << memory_Hessian     << "[Mb], "
        << "diagonalize: " << memory_eigenvector << "[Mb], "
        << "vectors: "     << _3N*_rfM*sizeof( double ) / 1000000 << "[Mb]."
        << endl;
      }
    }

  }

  void NormalModeDiagonalize::initialize( ProtoMolApp *app )
  {
    MTSIntegrator::initialize( app );

    // Setup adaptive timestep data from checkpoint
    app->eigenInfo.myOrigCEigval = origCEigVal;
    app->eigenInfo.myOrigTimestep = origTimestep;

    //test topmost integrator
    if ( top() != this ) {
      report << error << "NormalModeDiagonalize not top integrator." << endr;
    }

    myNextNormalMode  = dynamic_cast<NormalModeUtilities*>( myNextIntegrator );

    //Using complement of next integrator, so copy
    firstMode = myNextNormalMode->firstMode;
    numMode   = myNextNormalMode->numMode;

    //NM initialization, but fix dof, as set by next integrator, use complementry forces and dont generate noise
    NormalModeUtilities::initialize( ( int )app->positions.size(), app, myForces, COMPLIMENT_FORCES );
    _rfM = myNextNormalMode->_rfM;

    myLastNormalMode  = dynamic_cast<NormalModeUtilities*>( bottom() );

    //do first force calculation, and remove non sub-space part
    app->energies.clear(); //Need this or initial error, due to inner integrator energy?
    initializeForces();

    //Initialize Hessian array, OR assign hessian array for residues.
    if ( fullDiag ) {
      rHsn.initialData( _3N );
    } else {

      //automatically generate parameters?
      if(autoParmeters){
        residuesPerBlock = (int)pow((double)_N,0.6) / 15;
        blockVectorCols = 10 + (int)sqrt((float)residuesPerBlock);
        blockCutoffDistance = rHsn.cutOff;
        report << debug(1) << "[NormalModeDiagonalize::initialize] Auto parameters: residuesPerBlock " << residuesPerBlock <<
                        ", blockVectorCols " << blockVectorCols <<
                        ", blockCutoffDistance " << blockCutoffDistance << "." << endr;
      }

      //assign hessian array for residues, and clear.
      rHsn.initialResidueData( app->topology, residuesPerBlock, ( blockCutoffDistance == 0.0 ) );
    }

    // Check if array is already assigned
    if ( app->eigenInfo.myEigenvectors ) { 
      firstDiag = false;
      validMaxEigv = true;
    } else {
      firstDiag = true;
      validMaxEigv = false;

      //Calculate array size to be created
      app->eigenInfo.myEigenvectorLength = _N;
      app->eigenInfo.myNumEigenvectors = ( fullDiag == true ) ? _3N : _rfM;
      if(!app->eigenInfo.initializeEigenvectors())
          report << error << "Eigenvector array allocation error." << endr;

    }

    //flag used eigs
    app->eigenInfo.myNumUsedEigenvectors = _rfM;

    //Initialize BlockHessianDiagonalize, pass BlockHessian if Blocks (not full diag)
    if ( fullDiag ) {
      blockDiag.initialize( _3N );
    } else {
      blockDiag.initialize( &rHsn, _3N, (StandardIntegrator *)this );
    }

    //Diagnostics
    memory_Hessian = memory_eigenvector = 0;

    //setup rediag counter in case valid
    nextRediag = app->currentStep;//( int )( app->topology->time / getTimestep() ); //rediag first time

    //save positions where diagonalized for checkpoint save (assume I.C. if file)
    if(!checkpointUpdate) {
      diagAt = app->positions;
    } else {
      firstDiag = true;
    }

    newDiag = false;

    //timers/counters for diagnostics
    hessianCounter = rediagCounter = rediagUpdateCounter = 0;

  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  long NormalModeDiagonalize::run( const long numTimesteps ) {
    if ( numTimesteps < 1 ) return 0;
    
    //Current step at start
    int currentStepNum = app->currentStep; 

    //main loop
    app->energies.clear();
    for ( int i = 0;i < numTimesteps; ) {

      //Diagonalization if repetitive, first for forced
      if ( ( rediagCount && currentStepNum >= nextRediag ) || firstDiag || app->eigenInfo.reDiagonalize) {

        nextRediag += rediagCount;

        newDiag = true;

        report << debug(2) << "[NormalModeDiagonalize::run] Finding diagonalized Hessian." << endr;

        if(!(checkpointUpdate && firstDiag)) {
          //save positions where diagonalized for checkpoint save
          diagAt = app->positions;

          //remove last random perturbation?
          if ( removeRand ) {
            diagAt.intoSubtract( myLastNormalMode->gaussRandCoord1 );
          }
        }
        
        //save positions prior to diagonalze, if coarse then we need to use the 
        //  actual positions in case numeric S, but will be consistent for full (RJ)
        Vector3DBlock current_pos = app->positions;
        app->positions = diagAt;

        //~~~~if parallel only do diagonalization if master node~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        Real max_eig;
        
#ifdef HAVE_MPI
        if(Parallel::isParallel()){
          if(Parallel::iAmMaster()){
            //I am the MASTER if I get here
            
            //do actual diagonalization
            max_eig = doDiagonalization();

            //set new max eigenvalue in C
            app->eigenInfo.myNewCEigval = fabs(blockDiag.eigVal[_rfM]); //safe as eigval set t length sz=_3N >= _rfM

          }
          
          cout << "Diagonalized at master" << endl;
          
          //broadcast *Q from master
          const unsigned int numreals = app->eigenInfo.myEigenvectorLength * app->eigenInfo.myNumEigenvectors * 3;
          Parallel::bcastSlaves(*Q,*Q+numreals);
          
          //broadcast *eigValP from master
          Parallel::bcastSlaves(&max_eig,&max_eig + 1);
          
          //broadcast myNewCEigval from master
          Parallel::bcastSlaves(&app->eigenInfo.myNewCEigval,&app->eigenInfo.myNewCEigval + 1);
          
        }else{
#endif
          //do actual diagonalization
          max_eig = doDiagonalization();
          
          //set new max eigenvalue in C
          app->eigenInfo.myNewCEigval = fabs(blockDiag.eigVal[_rfM]); //safe as eigval set t length sz=_3N >= _rfM
#ifdef HAVE_MPI
        }
#endif
        
        //flag update to eigenvectors
        *eigVecChangedP = true;

        //set flags if firstDiag
        if ( firstDiag ) {
          numEigvectsu = ( fullDiag == true ) ? _3N : _rfM;
          *eigValP = max_eig;

          //first max eigenvalue in C, save original timestep for adaptive use
          if(!checkpointUpdate){
            app->eigenInfo.myOrigCEigval = app->eigenInfo.myNewCEigval;
            app->eigenInfo.myOrigTimestep = bottom()->getTimestep();
          }

          validMaxEigv = true;
          firstDiag = false;
        }

        //revert positions after diagonalization
        app->positions = current_pos;

        //post diag minimize?
        if(postDiagonalizeMinimize){
          Real lastLambda; int forceCalc = 0; //diagnostic/effective gamma
          
          app->eigenInfo.havePositionsChanged = true;

          //do minimization with local forces, max loop maxMinSteps, set subSpace minimization true
          int itrs = minimizer(minLim, maxMinSteps, true, false, true, &forceCalc, &lastLambda, &app->energies, &app->positions, app->topology, false, 0.0);
          report << debug(2) << "[NormalModeDiagonalize::run] iterations = "<< itrs << " force calcs = " << forceCalc << endr;

        }

        //adaptive timestep?
        const double newCEig = app->eigenInfo.myNewCEigval;
        const double oldCEig = app->eigenInfo.myOrigCEigval;
        const double baseTimestep = app->eigenInfo.myOrigTimestep;
        
        if(adaptiveTimestep && baseTimestep > 0 && newCEig > 0 && newCEig != oldCEig){
          
          const double tRatio = sqrt(oldCEig / newCEig);
          
          const double oldTimestep = bottom()->getTimestep();
    
          if (baseTimestep * tRatio <= baseTimestep) {
                ((STSIntegrator*)bottom())->setTimestep(baseTimestep * tRatio);
                
                report << debug(1) << "Adaptive time-step change, base " << baseTimestep << 
                                  ", new " << baseTimestep * tRatio << ", old " << oldTimestep << "." << endr;
          }
          
        }

        //sift current velocities/forces
        myNextNormalMode->subSpaceSift( &app->velocities, myForces );

        //clear re-diag flag
        app->eigenInfo.reDiagonalize = false;

      }

      //run integrator
      int stepsToRun = std::min((long)(nextRediag - currentStepNum), numTimesteps - i);
      
      const long completed = myNextIntegrator->run( stepsToRun );
      i += completed;
      currentStepNum += completed;

      //remove diagonalization flags after inner integrator call
      newDiag = false;
    }

    return numTimesteps;
  }

  //actual diagonalization code
  Real NormalModeDiagonalize::doDiagonalization(){
    
    //Diagonalize
    if ( fullDiag ) {
      //****Full method**********************************************************************//
      // Uses BLAS/LAPACK to do 'brute force' diagonalization                                //
      //*************************************************************************************//
      report << debug(2) << "Start diagonalization." << endr;
      
      //Find Hessians
      blockDiag.hessianTime.start(); //time Hessian
      rHsn.clear();
      rHsn.evaluate( &app->positions, app->topology, true ); //mass re-weighted hessian
      report << debug(2) << "Hessian found." << endr;
      
      //stop timer
      blockDiag.hessianTime.stop();
      hessianCounter++;
      
      //Diagonalize
      blockDiag.rediagTime.start();
      int numeFound;
      int info = blockDiag.diagHessian( *Q , blockDiag.eigVal, rHsn.hessM, _3N, numeFound );
      
      if ( info ) {
        report << error << "Full diagonalization failed." << endr;
      }
      
      //find number of -ve eigs
      int ii;
      for ( ii = 0;ii < _3N - 3;ii++ ) {
        if ( blockDiag.eigVal[ii+3] > 0 ) {
          break;
        }
      }
      
      report << debug( 1 ) << "[NormalModeDiagonalize::run] Full diagonalize. No. negative eigenvales = " << ii << endr;
      
      for ( int i = 0;i < _3N;i++ ) {
        blockDiag.eigIndx[i] = i;
      }
      
      blockDiag.absSort( *Q , blockDiag.eigVal, blockDiag.eigIndx, _3N );
      
      blockDiag.rediagTime.stop();
      rediagCounter++;
      
      //return max eig
      return blockDiag.eigVal[_3N-1];
      
    } else {
      //****Coarse method**************************************************************************//
      // Process:  Finds isolated 'minimized' block (of residues) Hessians   [evaluateResidues]    //
      //           Diagonalizes blocks to form block eigenvectors B          [findCoarseBlockEigs] //
      //           Finds actual Hessian H (but coarse grained) then S=B^THB  [innerHessian]        //
      //           Diagonalizes S to get eigenvectors Q, then approximate                          //
      //           eigenvectors are the first 'm' columns of BQ.                                   //
      //*******************************************************************************************//
      report << debug(2) << "Start coarse diagonalization." << endr;
      
      //find eienstuff
      Real max_eigenvalue = blockDiag.findEigenvectors( &app->positions, app->topology,
                                                       *Q , _3N, _rfM,
                                                       blockCutoffDistance, eigenValueThresh, blockVectorCols,
                                                       geometricfdof, numerichessians);
      
      //Stats/diagnostics
      rediagCounter++; hessianCounter++;
      memory_Hessian = ( rHsn.memory_base + rHsn.memory_blocks ) * sizeof( Real ) / 1000000;
      memory_eigenvector = blockDiag.memory_footprint * sizeof( Real ) / 1000000;
      
      //set new max eigenvalue in C
      app->eigenInfo.myNewCEigval = fabs(blockDiag.eigVal[_rfM]); //safe as eigval set t length sz=_3N >= _rfM
      
      //flag update to eigenvectors
      *eigVecChangedP = true;
            
      report << debug(2) << "Coarse diagonalization complete. Maximum eigenvalue = " << max_eigenvalue << "." << endr;
      
      //return max eig
      return max_eigenvalue;
    }
    
  }
  
  //********************************************************************************************************************************************

  //*************************************************************************************
  //****Output int paramiters************************************************************
  //*************************************************************************************

  void NormalModeDiagonalize::getParameters( vector<Parameter>& parameters ) const
  {
    MTSIntegrator::getParameters( parameters );

    parameters.push_back( Parameter( "reDiagFrequency",
                                     Value( rediagCount, ConstraintValueType::NotNegative() ),
                                     0,
                                     Text( "Frequency of re-diagonalization (steps)."        ) ) );

    parameters.push_back( Parameter( "fullDiag",
                                     Value( fullDiag,            ConstraintValueType::NoConstraints() ),
                                     false,
                                     Text( "Full diagonalization?"                           ) ) );

    parameters.push_back( Parameter( "removeRand",
                                     Value( removeRand,          ConstraintValueType::NoConstraints() ),
                                     false,
                                     Text( "Remove last random perturbation?"                ) ) );

    parameters.push_back( Parameter( "rediagHysteresis",
                                     Value( rediagHysteresis,    ConstraintValueType::NotNegative()   ),
                                     0.0,
                                     Text( "Re-diagonalization hysteresis."                  ) ) );

    parameters.push_back( Parameter( "eigenValueThresh",
                                     Value( eigenValueThresh,    ConstraintValueType::NotNegative()   ),
                                     5.0,
                                     Text( "'Inner' eigenvalue inclusion threshold."         ) ) );

    parameters.push_back( Parameter( "blockVectorCols",
                                     Value( blockVectorCols,     ConstraintValueType::NotNegative()   ),
                                     0,
                                     Text( "Target number of block eigenvector columns."     ) ) );

    parameters.push_back( Parameter( "residuesPerBlock",
                                     Value( residuesPerBlock,    ConstraintValueType::NotNegative()   ),
                                     1,
                                     Text( "Residues per block."                             ) ) );

    parameters.push_back( Parameter( "blockCutoffDistance",
                                     Value( blockCutoffDistance, ConstraintValueType::NotNegative()   ),
                                     10,
                                     Text( "Block cutoff distance for electrostatic forces." ) ) );

    parameters.push_back( Parameter( "autoParameters",
                                     Value(autoParmeters, ConstraintValueType::NoConstraints()   ),
                                     false, 
                                     Text("Automatically generate diagonalization parameters.") ) );

    parameters.push_back( Parameter( "adaptiveTimestep",
                                    Value(adaptiveTimestep, ConstraintValueType::NoConstraints()   ),
                                    false, 
                                    Text("Adapt time-step to latest diagonalization eigenvalues.") ) );

    parameters.push_back( Parameter( "postDiagonalizeMinimize",
                                    Value(postDiagonalizeMinimize, ConstraintValueType::NoConstraints()   ),
                                    false, 
                                    Text("Minimize after diagonalization.") ) );

    parameters.push_back( Parameter( "minimlim",
                                    Value(minLim,ConstraintValueType::NotNegative() ),
                                    0.1,
                                    Text("Minimizer target PE difference kcal mole^{-1}") ) );

    parameters.push_back( Parameter("maxminsteps",
                                    Value(maxMinSteps,ConstraintValueType::Positive()),
                                    100,
                                    Text("maximum number of minimizer steps.")));
      
    parameters.push_back( Parameter("geometricfdof",
                                    Value(geometricfdof, ConstraintValueType::NoConstraints()),
                                    false, Text("Calculate fixed degrees of freedom geometrically.")));
      
    parameters.push_back( Parameter("numericHessians",
                                    Value(numerichessians, ConstraintValueType::NoConstraints()),
                                    false, Text("Calculate Hessians numerically.")));


    
      }


  MTSIntegrator* NormalModeDiagonalize::doMake( const vector<Value>& values, ForceGroup* fg, StandardIntegrator *nextIntegrator ) const
  {
    return new NormalModeDiagonalize( values[0], values[1], values[2],
                                      values[3], values[4], values[5],
                                      values[6], values[7], values[8], 
                                      values[9], values[10], values[11],
                                      values[12], values[13],values[14], values[15], 
                                      fg, nextIntegrator               );
  }

  //*************************************************************************************
  //****Minimizers virtual force calculation*********************************************
  //*************************************************************************************
  void NormalModeDiagonalize::utilityCalculateForces()
  {
    app->energies.clear();

    calculateForces();
  }

  //*************************************************************************************
  //****Checkpointing********************************************************************
  //*************************************************************************************
  void NormalModeDiagonalize::streamRead( std::istream& inStream ) {

    inStream >> diagAt;
 
    //adaptive timestep
    inStream >> origCEigVal;    
    inStream >> origTimestep;
    
    checkpointUpdate = true;
      
  }
  
  void NormalModeDiagonalize::streamWrite( std::ostream& outStream ) const {    
    outStream.precision(15);
    outStream << diagAt;
    
    //adaptive timestep
    outStream << std::endl << app->eigenInfo.myOrigCEigval;
    outStream << " " << app->eigenInfo.myOrigTimestep;
    
  }


}
