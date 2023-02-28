/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Patrick Nathan
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* barn3Dturb.cpp:
 * This example examines a turbulent flow around a barn. At the
 * inlet, either a block profile or a power law wind profile is imposed as
 * a Dirchlet velocity boundary condition, whereas at the outlet a
 * Dirichlet pressure condition is set by p=0 (i.e. rho=1).
 * At both sides and on the top of the domain, a slip boundary condition is
 * imposed. On the ground and the walls of the barn, the bouzidi or bounce back
 * boundary
 * (no-slip) is imposed.
 *
 * Flow in x-direction with velocity=0 at height y=0.
 *
 * This example is adapted from the turbulence/nozzle3d.cpp example.
 */

#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used
#include "olb3D.hh"     // Include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;

// Use Bouzidi boundary on the ground and on the walls of the barn
// instead of Bounce Back
#define USE_BOUZIDI

// Choose your turbulent model of choice
//#define NoTurbulence // use BGKdynamics
#define Smagorinsky    //default turbulence model
//#define RLB
//#define ConsitentStrainSmagorinsky
//#define ShearSmagorinsky
//#define Krause

#ifdef ShearSmagorinsky
typedef D3Q19<AV_SHEAR> DESCRIPTOR;
#else
typedef D3Q19<> DESCRIPTOR;
#endif

// Parameters for the simulation setup
const int N = 5;                 // resolution of the model, for RLB N>=5,
                                 // others N>=2, but N>=5 recommended
const int M = 1;                 // time discretization refinement
const int inflowProfileMode = 1; // block profile (mode=0),
                                 // power profile (mode=1)
const T maxPhysT = 200.;         // max. simulation time in s, SI unit
const double Cs = 0.15;          // coefficient for Smagorinsky model
const T Vmax = 13.;              // maximal inlet velocity


// Define velocity profile
template <typename T, typename _DESCRIPTOR>
class TurbulentVelocity3D : public AnalyticalF3D<T,T> {

protected:
  // block profile (mode=0), power profile (mode=1)
  int _mode;
  T y_ref;
  T u_ref;
  T perturb;

public:
  TurbulentVelocity3D( UnitConverter<T,_DESCRIPTOR> const& converter, int mode=0 ) : AnalyticalF3D<T,T>( 3 )
  {
    _mode = mode;
    y_ref = 1.;
    u_ref = converter.getCharLatticeVelocity();
    perturb = 0.15;

    this->getName() = "turbulentVelocity3d";
  };

  bool operator()( T output[], const BaseType<T> input[] ) override
  {
    T y = input[1];
    // block profile inititalization
    T u_calc = u_ref;
    // power law profile inititalization
    if ( _mode==1 ) {
      T alpha  = 0.16;

      u_calc = u_ref * pow(y/y_ref, alpha);
    }

    T a = -1., b = 1.;
    T nRandom_0 = rand()/( T )RAND_MAX*( b-a ) + a;
    T nRandom_1 = rand()/( T )RAND_MAX*( b-a ) + a;
    T nRandom_2 = rand()/( T )RAND_MAX*( b-a ) + a;

    output[0] = u_calc + perturb*u_calc*nRandom_0;
    output[1] = perturb*u_calc*nRandom_1;
    output[2] = perturb*u_calc*nRandom_2;
    return true;
  };
};


void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // Material numbers:
  //  - 0:out of the domain
  //  - 1:fluid
  //  - 2:slip BC
  //  - 3:velocity inflow BC
  //  - 4:pressure outflow BC
  //  - 5:no-slip BC
  double lattice_length = converter.getConversionFactorLength();
  Vector<T,3> rMin = superGeometry.getStatistics().getMinPhysR( 0 );
  Vector<T,3> rMax = superGeometry.getStatistics().getMaxPhysR( 0 );
  clout<<" rMin: "<<rMin[0]<<" "<<rMin[1]<<" "<<rMin[2]<<std::endl;
  clout<<" rMax: "<<rMax[0]<<" "<<rMax[1]<<" "<<rMax[2]<<std::endl;

  // Set material number for boundaries (2)...
  superGeometry.rename( 0,2,indicator );
  // ...and inner fluid domain (1)
  superGeometry.rename( 2,1,stlReader );
  superGeometry.clean();

  // Set material number for inflow (3)
  Vector<T,3> origin = rMin + 0.5*lattice_length;
  origin[0] = rMin[0] - 0.5*lattice_length;
  Vector<T,3> extend = rMax - rMin - lattice_length;
  extend[0] = lattice_length;
  clout<<" (inflow) origin: "<<origin[0]<<" "<<origin[1]<<" "<<origin[2]<<std::endl;
  clout<<" (inflow) extend: "<<extend[0]<<" "<<extend[1]<<" "<<extend[2]<<std::endl;
  IndicatorCuboid3D<T> inflow( extend,origin );
  superGeometry.rename( 2,3,inflow );

  // Set material number for outflow (4)
  origin[0] = rMax[0] - 0.5*lattice_length;
  clout<<" (outflow) origin: "<<origin[0]<<" "<<origin[1]<<" "<<origin[2]<<std::endl;
  clout<<" (outflow) extend: "<<extend[0]<<" "<<extend[1]<<" "<<extend[2]<<std::endl;
  IndicatorCuboid3D<T> outflow( extend,origin );
  superGeometry.rename( 2,4,outflow );

  // Set material number for walls (5)
  // ground:
  origin = rMin - 0.5*lattice_length;
  extend = rMax - origin + 0.5*lattice_length;
  extend[1] = lattice_length;
  clout<<" (ground) origin: "<<origin[0]<<" "<<origin[1]<<" "<<origin[2]<<std::endl;
  clout<<" (ground) extend: "<<extend[0]<<" "<<extend[1]<<" "<<extend[2]<<std::endl;
  IndicatorCuboid3D<T> ground( extend,origin );
  superGeometry.rename( 2,5,ground );
  // barn walls:
  origin[2] = rMin[2] + lattice_length;
  extend[1] = (rMax[1] - origin[1])/2;
  extend[2] = rMax[2] - origin[2] - lattice_length;
  clout<<" (barn) origin: "<<origin[0]<<" "<<origin[1]<<" "<<origin[2]<<std::endl;
  clout<<" (barn) extend: "<<extend[0]<<" "<<extend[1]<<" "<<extend[2]<<std::endl;
  IndicatorCuboid3D<T> barn( extend,origin );
  superGeometry.rename( 2,5,barn );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     STLreader<T>& stlReader,
                     SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics (fluid)
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  sLattice.defineDynamics( superGeometry.getMaterialIndicator({1, 3, 4}), &bulkDynamics );

  setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
  setInterpolatedPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 4);

  // Material=2 -->symmetry (sides and top)
  sLattice.defineDynamics( superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>() );
  setSlipBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 2);

  // Material=5 -->bouzidi or bounce back (walls)
#ifdef USE_BOUZIDI
  sLattice.defineDynamics( superGeometry, 5, &instances::getNoDynamics<T, DESCRIPTOR>() );
   setBouzidiZeroVelocityBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 5, stlReader);
#else
  sLattice.defineDynamics( superGeometry, 5, &instances::getBounceBack<T, DESCRIPTOR>() );
#endif

  clout << "Prepare Lattice ... OK" << std::endl;
}


void setBoundandInitConditions( UnitConverter<T,DESCRIPTOR> const&converter,
                        SuperLattice3D<T,DESCRIPTOR>& lattice, SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"setBoundandInitConditions" );
  clout << "Set Boundary Conditions ..." << std::endl;

  AnalyticalConst3D<T,T> rhoF( 1 );
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> uF( velocity );

  // Seeding of random fluctuations and definition of the velocity field
  srand( time( nullptr ) );
  TurbulentVelocity3D<T,DESCRIPTOR> uSol( converter, inflowProfileMode );

  // Initialize all values of distribution functions to their local equilibrium
  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3, 4});

  lattice.defineRhoU( bulkIndicator, rhoF, uSol );

  lattice.iniEquilibrium( bulkIndicator, rhoF, uSol );

  // Make the lattice ready for simulation
  lattice.initialize();

  clout << "Set Boundary Conditions ... OK" << std::endl;
}


void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                        SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // Seeding of random fluctuations and definition of the velocity field
  srand( time( nullptr ) );
  TurbulentVelocity3D<T,DESCRIPTOR> uSol( converter, inflowProfileMode );

  sLattice.defineU( superGeometry, 3, uSol);
}


void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer )
{

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "barn3s_turb" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes the vtk files
  if ( iT%converter.getLatticeTime( maxPhysT/100. )==0 ) {
    // Create the data-reading functors...
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
    BlockReduction3D2D<T> planeReduction( normVel,Vector<T,3>({2., 0.5, 1}),
                                                  Vector<T,3>({0, 0, 1}));
    // write output as JPEG
    heatmap::write(planeReduction, iT);
  }

  // Writes output on the console
  if ( iT%converter.getLatticeTime( maxPhysT/200. )==0 ) {
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );
  }
}


int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  // clout.setMultiOutput(true);

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
     int {N},       // resolution: number of voxels per charPhysL
    (T)   0.500018, // latticeRelaxationTime: relaxation time, have to be
                    // greater than 0.5!
    (T)   0.2,      // charPhysLength: reference length of simulation geometry
    (T)   {Vmax},   // charPhysVelocity: maximal/highest expected velocity
                    // during simulation in __m / s__
    (T)   0.000015, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.204     // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("barn_turb");


  // === 2nd Step: Prepare Geometry ===

  STLreader<T> stlReader( "StallDomainFinal.stl", converter.getConversionFactorLength(), 0.001, 0, true);
  stlReader.print();
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getConversionFactorLength() );
  clout << "ConversionFactorLength " << converter.getConversionFactorLength() << std::endl;

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getConversionFactorLength(), noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );
  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );


  // === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  Dynamics<T, DESCRIPTOR>* bulkDynamics;
  const T omega = converter.getLatticeRelaxationFrequency();
#if defined(NoTurbulence)
  bulkDynamics = new BGKdynamics<T, DESCRIPTOR>( omega,
    instances::getBulkMomenta<T, DESCRIPTOR>() );
#elif defined(Smagorinsky)
  bulkDynamics = new SmagorinskyBGKdynamics<T, DESCRIPTOR>( omega,
    instances::getBulkMomenta<T, DESCRIPTOR>(), double {Cs});
#elif defined(RLB)
  bulkDynamics = new RLBdynamics<T, DESCRIPTOR>( omega,
    instances::getBulkMomenta<T, DESCRIPTOR>() );
#elif defined(ShearSmagorinsky)
  bulkDynamics = new ShearSmagorinskyBGKdynamics<T, DESCRIPTOR>( omega,
    instances::getBulkMomenta<T, DESCRIPTOR>(), double {Cs});
#elif defined(Krause)
  bulkDynamics = new KrauseBGKdynamics<T, DESCRIPTOR>( omega,
    instances::getBulkMomenta<T, DESCRIPTOR>(), double {Cs});
#else //ConsitentStrainSmagorinsky
  bulkDynamics = new ConStrainSmagorinskyBGKdynamics<T, DESCRIPTOR>( omega,
    instances::getBulkMomenta<T, DESCRIPTOR>(), double {Cs});
#endif

  prepareLattice(sLattice, converter, *bulkDynamics, stlReader, superGeometry);


  // === 4th Step: Set Boundary and Initial Conditions ===
  setBoundandInitConditions( converter, sLattice, superGeometry);


  // === 5th Step: Main Loop with Timer ===

  Timer<T> timer( converter.getLatticeTime( maxPhysT ),
    superGeometry.getStatistics().getNvoxel() );
  timer.start();

//////////////////////////
//  SuperVTMwriter3D<T> vtmWriter( "barn3d_turb" );
//  SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
//    vtmWriter.write( geometry );
/////////////////////////

  for ( std::size_t iT = 0; iT <= converter.getLatticeTime( maxPhysT ); ++iT ) {
//    // Apply filter (Approximate Deconvolution Method) simulating turbulence
//#ifdef ADM
//    SuperLatticeADM3D<T, DESCRIPTOR> admF( sLattice, 0.01, 2 );
//    admF.execute( superGeometry, 1 );
//#endif

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry );


    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();


    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer );

  }
  timer.stop();
  timer.printSummary();
  delete bulkDynamics;
}
