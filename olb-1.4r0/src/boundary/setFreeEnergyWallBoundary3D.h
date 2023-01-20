/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

///This file contains the Free Energy Wall Boundary
///This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_FREE_ENERGY_WALL_BOUNDARY_3D_H
#define SET_FREE_ENERGY_WALL_BOUNDARY_3D_H

#include <vector>
#include "utilities/functorPtr.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "geometry/superGeometry3D.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "core/superLattice3D.h"
#include "functors/lattice/indicator/superIndicatorBaseF3D.h"
#include "dynamics/dynamics.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "momentaOnBoundaries3D.h"
#include "io/ostreamManager.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/freeEnergyDynamics.h"
#include "boundaryPostProcessors3D.h"
#include "wallFunctionBoundaryPostProcessors3D.h"


namespace olb {
/// Implementation of a wetting boundary condition for the ternary free energy model, consisting of a BounceBack
/// dynamics and an FreeEnergyWall PostProcessor.
/// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
/// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
/// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly

////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry, int material,
                               T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber);

///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                               T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber);




/// Implementation of a wetting boundary condition for the ternary free energy model, consisting of a BounceBack
/// dynamics and an FreeEnergyWall PostProcessor.
/// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
/// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa3_ - Parameter related to surface tension. [lattice units]
/// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h3_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly

///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry, int material,
                               T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2,T h3, int latticeNumber);

///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                               T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber);



////////// BlockLattice Domain  /////////////////////////////////////////

//set FreeEnergyWallBoundary on block domain.
//This function works for the setFreeEnergyWallBoundary with h1,h2,h3 Parameters and h1,h2 Parameters
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator,
                               T addend, int latticeNumber, bool includeOuterCells=false);


}

#endif
