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
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systemes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2021/09/02
 *
 * Header file for module SplineCorrectedPolygonalCalculus.cpp
 *
 * This file is part of the DGtal library.
 */
//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include "DGtal/base/ConstAlias.h"
#include "DGtal/base/Common.h"
#include "DGtal/shapes/SurfaceMesh.h"
#include "DGtal/math/linalg/EigenSupport.h"
#include "PolygonalCalculus.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// template class SplineCorrectedPolygonalCalculus
/**
 * Description of template class 'SplineCorrectedPolygonalCalculus' <p>
 * \brief Implements differential operators on polygonal surfaces from
 * @cite degoes2020discrete
 *
 * See @ref moduleSplineCorrectedPolygonalCalculus for details.
 *
 * @tparam TRealPoint a model of points R^3 (e.g. PointVector).
 * @tparam TRealVector a model of vectors in R^3 (e.g. PointVector).
 */
template <typename TRealPoint, typename TRealVector>
class SplineCorrectedPolygonalCalculus : public PolygonalCalculus<TRealPoint,TRealVector>
{
    // ----------------------- Standard services ------------------------------
public:

    ///Concept checking
    static const Dimension dimension = TRealPoint::dimension;
    BOOST_STATIC_ASSERT( ( dimension == 3 ) );

    ///parent type
    typedef PolygonalCalculus<TRealPoint,TRealVector> Calculus;

    ///Self type
    typedef SplineCorrectedPolygonalCalculus<TRealPoint, TRealVector> Self;

    ///Type of SurfaceMesh
    typedef SurfaceMesh<TRealPoint, TRealVector> MySurfaceMesh;
    ///Vertex type
    typedef typename MySurfaceMesh::Vertex Vertex;
    ///Face type
    typedef typename MySurfaceMesh::Face Face;
    ///Position type
    typedef typename MySurfaceMesh::RealPoint Real3dPoint;
    ///Real vector type
    typedef typename MySurfaceMesh::RealVector Real3dVector;


    ///Linear Algebra Backend from Eigen
    typedef EigenLinearAlgebraBackend LinAlg;
    ///Type of Vector
    typedef LinAlg::DenseVector Vector;
    ///Type of dense matrix
    typedef LinAlg::DenseMatrix DenseMatrix;
    ///Type of sparse matrix
    typedef LinAlg::SparseMatrix SparseMatrix;
    ///Type of sparse matrix triplet
    typedef LinAlg::Triplet Triplet;

    ///Type of a sparse matrix solver
    typedef LinAlg::SolverSimplicialLDLT Solver;

    typedef Eigen::Vector3d Vector3;

    SplineCorrectedPolygonalCalculus(
              const ConstAlias<MySurfaceMesh> surf,
              const std::function<Vector(Vertex)> &normal_embedder,
              bool globalInternalCacheEnabled = false):
        PolygonalCalculus<TRealPoint,TRealVector>(surf,normal_embedder,globalInternalCacheEnabled)
    {
    };

    SplineCorrectedPolygonalCalculus(
              const ConstAlias<MySurfaceMesh> surf,
              const std::function<Real3dPoint(Face,Vertex)> &pos_embedder,
              const std::function<Vector(Vertex)> &normal_embedder,
              bool globalInternalCacheEnabled = false):
        PolygonalCalculus<TRealPoint,TRealVector>(surf,pos_embedder,normal_embedder,globalInternalCacheEnabled)
    {
    };

    /**
   * Deleted default constructor.
   */
    SplineCorrectedPolygonalCalculus() = delete;

    /**
   * Destructor (default).
   */
    ~SplineCorrectedPolygonalCalculus() = default;

    /**
   * Deleted copy constructor.
   * @param other the object to clone.
   */
    SplineCorrectedPolygonalCalculus ( const SplineCorrectedPolygonalCalculus & other ) = delete;

    /**
   * Deleted move constructor.
   * @param other the object to move.
   */
    SplineCorrectedPolygonalCalculus ( SplineCorrectedPolygonalCalculus && other ) = delete;

    /**
   * Deleted copy assignment operator.
   * @param other the object to copy.
   * @return a reference on 'this'.
   */
    SplineCorrectedPolygonalCalculus & operator= ( const SplineCorrectedPolygonalCalculus & other ) = delete;

    /**
   * Deleted move assignment operator.
   * @param other the object to move.
   * @return a reference on 'this'.
   */
    SplineCorrectedPolygonalCalculus & operator= ( SplineCorrectedPolygonalCalculus && other ) = delete;

    // ----------------------- Inner Spline Classes --------------------------------------
public:
struct Spline{
    static constexpr int DEGREE = 3;
    DenseMatrix coeffs;

    Vector midpointForm;

    Spline() {
        midpointForm = Vector::Zero(4);
        midpointForm(1) = 0.5;
        midpointForm(2) = 1./3.;
        midpointForm(3) = 0.25;
    }

    Vector eval(double t) const {
        Vector T = Vector::Zero(4);
        double x = 1;
        for (int i = 0;i<4;i++){
            T(i) = x;
            x *= t;
        }
        return (coeffs*T).col(0);
    }

    Vector operator()(double t) const{ return eval(t);}

    TRealPoint evalAsPoint(double t) const {
        Vector X = eval(t);
        TRealPoint Xp;
        for (int i = 0;i<3;i++)
            Xp(i) = X(i);
        return Xp;
    }
    Vector evalTangent(double t) const {
        Vector T = Vector::Zero(4);
        T(3) = 3*t*t;
        T(2) = 2*t;
        T(1) = 1;
        return (coeffs*T).col(0).normalized();
    }
    TRealPoint evalTangentAsPoint(double t) const {
        Vector X = evalTangent(t);
        TRealPoint Xp;
        for (int i = 0;i<3;i++)
            Xp(i) = X(i);
        return Xp;
    }
    Vector evalNormal(double t) const {
        Vector T = Vector::Zero(4);
        T(3) = 6*t;
        T(2) = 2;
        return (coeffs*T).col(0).normalized();
    }
    TRealPoint evalNormalAsPoint(double t) const {
        Vector X = evalNormal(t);
        TRealPoint Xp;
        for (int i = 0;i<3;i++)
            Xp(i) = X(i);
        return Xp;
    }
    Vector averagePoint() const {
        return coeffs*midpointForm;
    }
};

struct SplineMaker {

    Eigen::ColPivHouseholderQR<Eigen::Matrix4d> normal_solver;
    Eigen::ColPivHouseholderQR<Eigen::Matrix4d> tangent_solver;
    SplineMaker(){
        Eigen::Matrix4d N;
        N << 1,0,0,0,
             1,1,1,1,
             0,0,2,0,
             0,0,2,6;
        normal_solver.compute(N);
        Eigen::Matrix4d T;
        T << 1,0,0,0,
             1,1,1,1,
             0,1,0,0,
             0,1,2,3;
        tangent_solver.compute(T);
    }
    Spline makeNormalSpline(
            const Vector& x1,
            const Vector& n1,
            const Vector& x2,
            const Vector& n2) const{
        Spline S;
        DenseMatrix B(4,3);
        for (int i =0;i<3;i++){
            B(0,i) = x1(i);
            B(1,i) = x2(i);
            B(2,i) = -n1(i);
            B(3,i) = -n2(i);
        }
        S.coeffs = normal_solver.solve(B).transpose();
        return S;
    }
    Spline makeTangentSpline(
            const Vector& x1,
            const Vector& t1,
            const Vector& x2,
            const Vector& t2) const{
        Spline S;
        DenseMatrix B(4,3);
        for (int i =0;i<3;i++){
            B(0,i) = x1(i);
            B(1,i) = x2(i);
            B(2,i) = t1(i);
            B(3,i) = t2(i);
        }
        S.coeffs = tangent_solver.solve(B).transpose();
        return S;
    }
};

    SplineMaker splineMaker;

public:
    // ---------------------- Redefined Operators ---------------------------

    static Eigen::Vector3d toVec3(const Vector& x){
        return Eigen::Vector3d(x(0),x(1),x(2));
    }
    static Eigen::Vector3d toVec3(const Real3dPoint& x){
        return Eigen::Vector3d(x(0),x(1),x(2));
    }

    /// Polygonal (corrected) vector area.
    /// @param f the face
    /// @return a vector
    Vector vectorArea(const Face f) const override
    {
        auto vertices = this->mySurfaceMesh->incidentVertices(f);
        auto nf = vertices.size();
        Vector af = Vector::Zero(3);
        static const double lambda = 1./60.;
        for (auto v = 0u;v<nf;v++){
            auto i = vertices[v];
            Vector3 xi = toVec3(this->myEmbedder(f,i));
            Vector3 ni = toVec3(this->myVertexNormalEmbedder(i));
            auto j = vertices[(v+1)%nf];
            Vector3 xj = toVec3(this->myEmbedder(f,j));
            auto nj = toVec3(this->myVertexNormalEmbedder(j));
            Spline S = splineMaker.makeNormalSpline(
                        xi,
                        ni,
                        xj,
                        nj
                        );
            Vector3 T0 = toVec3(S.evalTangent(0));
            Vector3 T1 = toVec3(S.evalTangent(1));
            Vector3 TH = toVec3(S.evalTangent(0.5));
            af += lambda*(
                        14*(xi.cross(T0) + xj.cross(T1))
                        +32*toVec3(S(0.5)).cross(TH)
                        + xi.cross(ni)
                        - xj.cross(nj)
                        );
        }
        return af;
    }

    DenseMatrix coGradient(const Face f) const override
    {
        auto vertices = this->mySurfaceMesh->incidentVertices(f);
        auto nf = vertices.size();
        DenseMatrix CGS = DenseMatrix::Zero(3,nf);
        static const double lambda = 1./6.;
        for (auto v = 0u;v<nf;v++){
            auto i = vertices[v];
            auto j = vertices[(v+1)%nf];
            Spline S = splineMaker.makeNormalSpline(
                        toVec3(this->myEmbedder(f,i)),
                        toVec3(this->myVertexNormalEmbedder(i)),
                        toVec3(this->myEmbedder(f,j)),
                        toVec3(this->myVertexNormalEmbedder(j))
                        );
            auto T0 = S.evalTangent(0);
            auto T1 = S.evalTangent(1);
            auto TH = S.evalTangent(0.5);
            CGS.col(v) += lambda*(T0 + 2*TH);
            CGS.col((v+1)%nf) += lambda*(T1 + 2*TH);
        }
        return CGS;
    }

    DenseMatrix B(const Face f) const override
    {
        auto vertices = this->mySurfaceMesh->incidentVertices(f);
        auto nf = vertices.size();
        DenseMatrix midpoints = DenseMatrix::Zero(nf,3);
        for (auto v = 0u;v<nf;v++){
            auto i = vertices[v];
            auto j = vertices[(v+1)%nf];
            Spline S = splineMaker.makeNormalSpline(
                        toVec3(this->myEmbedder(f,i)),
                        toVec3(this->myVertexNormalEmbedder(i)),
                        toVec3(this->myEmbedder(f,j)),
                        toVec3(this->myVertexNormalEmbedder(j))
                        );
            midpoints.block(v,0,1,3) = S(0.5).transpose();
            //midpoints.block(v,0,1,3) = S.getMidPoint().transpose();
        }
        return midpoints;
    }

    Vector projectOnVertexTangentPlane(const Vector& e,Vertex v){
        DenseMatrix T = this->Tv(v);
        return (T*T.transpose()*e).col(0);
    }
}; // end of class SplineCorrectedPolygonalCalculus

/**
 * Overloads 'operator<<' for displaying objects of class 'SplineCorrectedPolygonalCalculus'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'SplineCorrectedPolygonalCalculus' to write.
 * @return the output stream after the writing.
 */
template <typename TP, typename TV>
std::ostream&
operator<< ( std::ostream & out, const SplineCorrectedPolygonalCalculus<TP,TV> & object )
{
    object.selfDisplay( out );
    return out;
}

} // namespace DGtal
///////////////////////////////////////////////////////////////////////////////
