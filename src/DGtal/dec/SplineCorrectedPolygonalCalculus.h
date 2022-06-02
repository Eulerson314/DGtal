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
protected:
struct Spline{
    static constexpr int DEGREE = 3;
    DenseMatrix coeffs;

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
};

struct SplineMaker {

    Eigen::ColPivHouseholderQR<Eigen::Matrix4d> solver;
    SplineMaker(){
        Eigen::Matrix4d A = Eigen::Matrix4d::Zero(4,4);
        A(0,0) = 1;
        A.block<1,4>(1,0) = Eigen::MatrixXd::Ones(1,4);
        A(2,2) = 2;
        A(3,2) = 2; A(3,3) = 6;
        solver.compute(A);
    }
    Spline makeSpline(
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
        S.coeffs = solver.solve(B).transpose();
        return S;
    }
};

    SplineMaker splineMaker;

public:
    // ---------------------- Redefined Operators ---------------------------

    /// Polygonal (corrected) vector area.
    /// @param f the face
    /// @return a vector
    Vector vectorArea(const Face f) const override
    {
        auto vertices = mySurfaceMesh->incidentVertices(f);
        auto nf = vertices.size();
        Vector af = Vector::Zero(3);
        static const double lambda = 1./60.;
        for (auto v = 0u;v<nf;v++){
            auto i = vertices[v];
            Vector3 xi = Calculus::toVec3(myEmbedder(f,i));
            Vector3 ni = toVec3(myVertexNormalEmbedder(i));
            auto j = vertices[(v+1)%nf];
            Vector3 xj = Calculus::toVec3(myEmbedder(f,j));
            auto nj = toVec3(myVertexNormalEmbedder(j));
            Spline S = splineMaker.makeSpline(
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
        auto vertices = mySurfaceMesh->incidentVertices(f);
        auto nf = vertices.size();
        DenseMatrix CGS = DenseMatrix::Zero(3,nf);
        static const double lambda = 1./6.;
        for (auto v = 0u;v<nf;v++){
            auto i = vertices[v];
            auto j = vertices[(v+1)%nf];
            Spline S = splineMaker.makeSpline(
                        Calculus::toVec3(myEmbedder(f,i)),
                        toVec3(myVertexNormalEmbedder(i)),
                        Calculus::toVec3(myEmbedder(f,j)),
                        toVec3(myVertexNormalEmbedder(j))
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
        auto vertices = mySurfaceMesh->incidentVertices(f);
        auto nf = vertices.size();
        DenseMatrix midpoints = DenseMatrix::Zero(nf,3);
        for (auto v = 0u;v<nf;v++){
            auto i = vertices[v];
            auto j = vertices[(v+1)%nf];
            Spline S = splineMaker.makeSpline(
                        Calculus::toVec3(myEmbedder(f,i)),
                        toVec3(myVertexNormalEmbedder(i)),
                        Calculus::toVec3(myEmbedder(f,j)),
                        toVec3(myVertexNormalEmbedder(j))
                        );
            midpoints.block(0,v,3,1) = S(0.5);
        }
        return midpoints;
    }

    DenseMatrix M(const Face f, const double lambda=1.0) const override
    {
        if (checkCache(M_,f))
            return myGlobalCache[M_][f];

        auto Gf=gradient(f);
        auto Pf=P(f);
        DenseMatrix op = faceArea(f) * Gf.transpose()*Gf + lambda * Pf.transpose()*Pf;

        setInCache(M_,f,op);
        return op;
    }



    // ----------------------- Interface --------------------------------------


    /// Update the embedding function.
    /// @param externalFunctor a new embedding functor (Face,Vertex)->RealPoint.
    void setEmbedder(const std::function<Real3dPoint(Face,Vertex)> &externalFunctor)
    {
        myEmbedder = externalFunctor;
    }

    // ----------------------- Cache mechanism --------------------------------------

    /// Generic method to compute all the per face DenseMatrices and store them in an
    /// indexed container.
    ///
    /// Usage example:
    /// @code
    ///auto opM = [&](const SplineCorrectedPolygonalCalculus<Mesh>::Face f){ return calculus.M(f);};
    ///auto cacheM = boxCalculus.getOperatorCacheMatrix(opM);
    ///...
    /////Now you have access to the cached values and mixed them with un-cached ones
    ///  Face f = ...;
    ///  auto res = cacheM[f] * calculus.D(f) * phi;
    /// ...
    ///@endcode
    ///
    /// @param perFaceOperator the per face operator
    /// @return an indexed container of all DenseMatrix operators (indexed per Face).
    std::vector<DenseMatrix> getOperatorCacheMatrix(const std::function<DenseMatrix(Face)> &perFaceOperator) const
    {
        std::vector<DenseMatrix> cache;
        for(auto f=0; f < mySurfaceMesh->nbFaces(); ++f)
            cache.push_back(perFaceOperator(f));
        return cache;
    }

    /// Generic method to compute all the per face Vector and store them in an
    /// indexed container.
    ///
    /// Usage example:
    /// @code
    ///auto opCentroid = [&](const SplineCorrectedPolygonalCalculus<Mesh>::Face f){ return calculus.centroid(f);};
    ///auto cacheCentroid = boxCalculus.getOperatorCacheVector(opCentroid);
    ///...
    /////Now you have access to the cached values and mixed them with un-cached ones
    ///  Face f = ...;
    ///  auto res = calculus.P(f) * cacheCentroid[f] ;
    /// ...
    ///@endcode
    ///
    /// @param perFaceVectorOperator the per face operator
    /// @return an indexed container of all Vector quantities (indexed per Face).
    std::vector<Vector> getOperatorCacheVector(const std::function<Vector(Face)> &perFaceVectorOperator) const
    {
        std::vector<Vector> cache;
        for(auto f=0; f < mySurfaceMesh->nbFaces(); ++f)
            cache.push_back(perFaceVectorOperator(f));
        return cache;
    }

    /// Enable the internal global cache for operators.
    ///
    void enableInternalGlobalCache()
    {
        myGlobalCacheEnabled = true;
    }

    /// Disable the internal global cache for operators.
    /// This method will also clean up the
    void disableInternalGlobalCache()
    {
        myGlobalCacheEnabled = false;
        myGlobalCache.clear();
    }

    // ----------------------- Common --------------------------------------
public:

    /// Update the internal cache structures
    /// (e.g. degree of each face).
    void init()
    {
        updateFaceDegree();
    }

    /// Helper to retrieve the degree of the face from the cache.
    /// @param f the face
    /// @return the number of vertices of the face.
    size_t faceDegree(Face f) const
    {
        return myFaceDegree[f];
    }

    /// @return the number of vertices of the underlying surface mesh.
    size_t nbVertices() const
    {
        return mySurfaceMesh->nbVertices();
    }

    /// @return the number of faces of the underlying surface mesh.
    size_t nbFaces() const
    {
        return mySurfaceMesh->nbFaces();
    }

    /// @returns the degree of the face f (number of vertices)
    /// @param f the face
    size_t degree(const Face f) const
    {
        return myFaceDegree[f];
    }

    /// @returns an pointer to the underlying SurfaceMash object.
    const MySurfaceMesh * getSurfaceMeshPtr() const
    {
        return mySurfaceMesh;
    }

    /**
   * Writes/Displays the object on an output stream.
   * @param out the output stream where the object is written.
   */
    void selfDisplay ( std::ostream & out ) const
    {
        out << "[SplineCorrectedPolygonalCalculus]: ";
        if (myGlobalCacheEnabled)
            out<< "internal cache enabled, ";
        else
            out<<"internal cache disabled, ";
        out <<"SurfaceMesh="<<*mySurfaceMesh;
    }

    /**
   * Checks the validity/consistency of the object.
   * @return 'true' if the object is valid, 'false' otherwise.
   */
    bool isValid() const
    {
        return true;
    }

    // ------------------------- Protected Datas ------------------------------
protected:

    ///Enum for operators in the internal cache strategy
    enum OPERATOR { X_, D_, E_, A_, COGRAD_, GRAD_, FLAT_, B_, SHARP_, P_, M_, DIVERGENCE_, CURL_, L_ };

    /// Update the face degree cache
    void updateFaceDegree()
    {
        myFaceDegree.resize(mySurfaceMesh->nbFaces());
        for(auto f = 0; f <  mySurfaceMesh->nbFaces(); ++f)
        {
            auto vertices = mySurfaceMesh->incidentVertices(f);
            auto nf = vertices.size();
            myFaceDegree[f] = nf;
        }
    }

    /// Check internal cache if enabled.
    /// @param key the operator name
    /// @param f the face
    /// @returns true if the operator "key" for the face f has been computed.
    bool checkCache(OPERATOR key, const Face f) const
    {
        if (myGlobalCacheEnabled)
            if (myGlobalCache[key].find(f) != myGlobalCache[key].end())
                return true;
        return false;
    }

    /// Set an operator in the internal cache.
    /// @param key the operator name
    /// @param f the face
    /// @param ope the operator to store
    void setInCache(OPERATOR key, const Face f,
            const DenseMatrix &ope) const
    {
        if (myGlobalCacheEnabled)
            myGlobalCache[key][f]  = ope;
    }


    // ------------------------- Internals ------------------------------------
protected:


    ///Underlying SurfaceMesh
    const MySurfaceMesh *mySurfaceMesh;

    ///Embedding function (face,vertex)->R^3 for the vertex position wrt. the face.
    std::function<Real3dPoint(Face, Vertex)> myEmbedder;

    ///Embedding function (vertex)->R^3 for the vertex normal.
    std::function<Vector(Vertex)> myVertexNormalEmbedder;

    ///Cache containing the face degree
    std::vector<size_t> myFaceDegree;

    ///Global cache
    bool myGlobalCacheEnabled;
    mutable std::array<std::unordered_map<Face,DenseMatrix>, 14> myGlobalCache;

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
