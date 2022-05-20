/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file
 * @author Baptiste GENEST (\c baptistegenest@gmail.com )
 * Laboratoire d'InfoRmatique en Image et Systemes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2022/05/20
 *
 * Header file for module VectorsInHeat.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(VectorsInHeat_RECURSES)
#error Recursive header files inclusion detected in VectorsInHeat.h
#else // defined(VectorsInHeat_RECURSES)
/** Prevents recursive inclusion of headers. */
#define VectorsInHeat_RECURSES

#if !defined VectorsInHeat_h
/** Prevents repeated inclusion of headers. */
#define VectorsInHeat_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/base/ConstAlias.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{
/////////////////////////////////////////////////////////////////////////////
// template class VectorInHeat
/**
 * Description of template class 'VectorsInHeat' <p>
 * \brief This class implements @cite Crane2013 on polygonal surfaces  (using @ref modulePolygonalCalculus).
 *
 * see @ref moduleVectorsInHeat for details and examples.
 *
 * @tparam a model of PolygonalCalculus.
 */
template <typename TPolygonalCalculus>
class VectorsInHeat
{
    // ----------------------- Standard services ------------------------------
public:

    typedef TPolygonalCalculus PolygonalCalculus;
    typedef typename PolygonalCalculus::SparseMatrix SparseMatrix;
    typedef typename PolygonalCalculus::DenseMatrix DenseMatrix;
    typedef typename PolygonalCalculus::Solver Solver;
    typedef typename PolygonalCalculus::Vector Vector;
    typedef typename PolygonalCalculus::Vertex Vertex;

    /**
     * Default constructor.
     */
    VectorsInHeat() = delete;

    /// Constructor from an existing polygonal calculus. T
    /// @param calculus a instance of PolygonalCalculus
    VectorsInHeat(ConstAlias<PolygonalCalculus> calculus): myCalculus(&calculus)
    {
        myIsInit=false;
    }

    /**
     * Destructor.
     */
    ~VectorsInHeat() = default;

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    VectorsInHeat ( const VectorsInHeat & other ) = delete;

    /**
     * Move constructor.
     * @param other the object to move.
     */
    VectorsInHeat ( VectorsInHeat && other ) = delete;

    /**
     * Copy assignment operator.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    VectorsInHeat & operator= ( const VectorsInHeat & other ) = delete;

    /**
     * Move assignment operator.
     * @param other the object to move.
     * @return a reference on 'this'.
     */
    VectorsInHeat & operator= ( VectorsInHeat && other ) = delete;



    // ----------------------- Interface --------------------------------------

    /// Initialize the solvers with @a dt as timestep for the
    /// heat diffusion
    /// @param dt timestep
    void init(double dt)
    {
        myIsInit=true;

        //As the LB is PSD, the identity matrix shouldn't be necessary. However, some solvers
        //may have issues with positive semi-definite matrix.
        SparseMatrix I(myCalculus->nbVertices(),myCalculus->nbVertices());
        I.setIdentity();
        SparseMatrix laplacian = myCalculus->globalLaplaceBeltrami()  + I*1e-6;

        SparseMatrix I2(2*myCalculus->nbVertices(),2*myCalculus->nbVertices());
        I2.setIdentity();
        SparseMatrix connectionLaplacian = myCalculus->globalConnectionLaplace() + I2*1e-6;

        SparseMatrix mass      = myCalculus->globalLumpedMassMatrix();
        SparseMatrix mass2      = myCalculus->doubledGlobalLumpedMassMatrix();
        SparseMatrix scalarHeatOpe   =  mass + dt*laplacian;
        SparseMatrix vectorHeatOpe   =  mass2 + dt*connectionLaplacian;

        //Prefactorizing
        myScalarHeatSolver.compute(scalarHeatOpe);
        myVectorHeatSolver.compute(vectorHeatOpe);

        //empty sources
        myVectorSource     	= Vector::Zero(2*myCalculus->nbVertices());
        myScalarSource		= Vector::Zero(myCalculus->nbVertices());
        myDiracSource		= Vector::Zero(myCalculus->nbVertices());
    }

    /** Adds a source vector (3D extrinsic) at a vertex @e aV
    * @param aV the Vertex
     **/
    void addSource(const Vertex aV,const Vector& ev)
    {
        ASSERT_MSG(aV < myCalculus->nbVertices(), "Vertex is not in the surface mesh vertex range");
        Vector v = myCalculus->flat(aV)*myCalculus->sharp(aV)*ev;
        myVectorSource( 2*aV ) = v(0);
        myVectorSource( 2*aV+1 ) = v(1);
        myScalarSource( aV ) = v.norm();
        myDiracSource( aV ) = 1;
    }

    /**
     * @returns the source point vector.
     **/
    Vector vectorSource() const
    {
        FATAL_ERROR_MSG(myIsInit, "init() method must be called first");
        return myVectorSource;
    }

    Vector vectorSourceAtVertex(const Vertex aV){
        FATAL_ERROR_MSG(myIsInit, "init() method must be called first");
        DenseMatrix Tv = myCalculus->Tv(aV);
        return Tv.col(0)*myVectorSource(2*aV) + Tv.col(1)*myVectorSource(2*aV+1);
    }


    /// Main computation of the Vectors In Heat
    /// @returns the estimated heat diffused vectors from the sources.
    Vector compute() const
    {
        FATAL_ERROR_MSG(myIsInit, "init() method must be called first");
        //Heat diffusion
        Vector vectorHeatDiffusion = myVectorHeatSolver.solve(myVectorSource);
        Vector scalarHeatDiffusion = myScalarHeatSolver.solve(myScalarSource);
        Vector diracHeatDiffusion = myScalarHeatSolver.solve(myDiracSource);
        auto surfmesh = myCalculus->getSurfaceMeshPtr();

        Vector result(2*surfmesh->nbVertices());

        for (auto v = 0;v<surfmesh->nbVertices();v++){
            Vector Y(2);
            Y(0) = vectorHeatDiffusion(2*v);
            Y(1) = vectorHeatDiffusion(2*v+1);
            Y = Y.normalized()*(scalarHeatDiffusion(v)/diracHeatDiffusion(v));
            result(2*v) = Y(0);
            result(2*v+1) = Y(1);
        }

        return result;
    }


    /// @return true if the calculus is valid.
    bool isValid() const
    {
        return myIsInit && myCalculus->isValid();
    }

    // ----------------------- Private --------------------------------------

private:

    ///The underlying PolygonalCalculus instance
    const PolygonalCalculus *myCalculus;

    ///Poisson solver
    Solver myPoissonSolver;

    ///Heat solver
    Solver myScalarHeatSolver;
    Solver myVectorHeatSolver;

    ///Source vector
    Vector myScalarSource;
    Vector myDiracSource;
    Vector myVectorSource;

    ///Validitate flag
    bool myIsInit;

}; // end of class VectorsInHeat
} // namespace DGtal

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined VectorsInHeat_h

#undef VectorsInHeat_RECURSES
#endif // else defined(VectorsInHeat_RECURSES)
