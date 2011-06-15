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

/**
 * @file ctopo-2-3d.cpp
 * @ingroup Examples
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/05/22
 *
 * An example file named ctopo-2-3d.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"

#include "DGtal/helpers/StdDefs.h"
#include <QtGui/qapplication.h>
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/helpers/Surfaces.h"

#include "DGtal/io-viewers/3dViewers/DGtalQGLViewer.h"
#include "DGtal/io-viewers/readers/VolReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/io-viewers/colormaps/GradientColorMap.h"

#include "ConfigExamples.h"

///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  trace.beginBlock ( "Example ctopo-2-3d" );
  // for 3D display with DGtalQGLViewer
  QApplication application(argc,argv);
  
  typedef ImageSelector < Z3i::Domain, int>::Type Image;
  std::string inputFilename = examplesPath + "samples/cat10.vol"; 
  Image image = VolReader<Image>::importVol(inputFilename);
  Z3i::DigitalSet set3d (image.domain());
  SetFromImage<Z3i::DigitalSet>::append<Image>(set3d, image, 0,255);
  DGtalQGLViewer viewer;  
  viewer.show(); 
  
  
  // Construct the Khalimsky space from the image domain
  Z3i::KSpace ks;
  bool space_ok = ks.init( image.domain().lowerBound(), image.domain().upperBound(), true );
  
  std::vector<Z3i::SCell> vectBdrySCell;
  std::vector<Z3i::SCell> vectBdrySCell2;
  std::set<Z3i::SCell> vectBdrySCellALL;
  SurfelAdjacency<3> SAdj( true );
  

  
  //Extract an initial boundary cell
  Z3i::SCell aCell = Surfaces<Z3i::KSpace>::findABel(ks, set3d);
  
  // Extracting all boundary surfels which are connected to the initial boundary Cell.
  Surfaces<Z3i::KSpace>::trackBoundary( vectBdrySCellALL,
					ks,SAdj, set3d, aCell );
    
  // Extract the bondary contour associated to the initial surfel in its first direction
  Surfaces<Z3i::KSpace>::track2DBoundary( vectBdrySCell,
 					  ks, *(ks.sDirs( aCell )),SAdj, set3d, aCell );
  
  // Extract the bondary contour associated to the initial surfel in its second direction
  Surfaces<Z3i::KSpace>::track2DBoundary( vectBdrySCell2,
 					  ks, *(++(ks.sDirs( aCell ))),SAdj, set3d, aCell );  
  
  
  // Displaying all the surfels in transparent mode
  viewer << SetMode3D((*(vectBdrySCellALL.begin())).styleName(), "Transparent");
  for( std::set<Z3i::SCell>::iterator it=vectBdrySCellALL.begin(); 
       it!= vectBdrySCellALL.end(); it++){
    viewer<< *it;
  } 
  
  // Displaying First surfels cut with gradient colors.;
  GradientColorMap<int> cmap_grad( 0, vectBdrySCell2.size() );
  cmap_grad.addColor( DGtalBoard::Color( 50, 50, 255 ) );
  cmap_grad.addColor( DGtalBoard::Color( 255, 0, 0 ) );
  cmap_grad.addColor( DGtalBoard::Color( 255, 255, 10 ) );
  
  // Need to avoid surfel superposition (the surfel size in increased)
  viewer << DGtalQGLViewer::shiftSurfelVisu; 
  viewer << SetMode3D((*(vectBdrySCell2.begin())).styleName(), "");
  viewer.setFillColor(QColor(180, 200, 25, 255));
  
  int d=0;
  for( std::vector<Z3i::SCell>::iterator it=vectBdrySCell2.begin(); 
       it!= vectBdrySCell2.end(); it++){
    DGtalBoard::Color col= cmap_grad(d);
    viewer.setFillColor(QColor(col.red(),col.green() ,col.blue(), 255));
    viewer<< *it;
    d++;
  }
  
  GradientColorMap<int> cmap_grad2( 0, vectBdrySCell.size() );
  cmap_grad2.addColor( DGtalBoard::Color( 50, 50, 255 ) );
  cmap_grad2.addColor( DGtalBoard::Color( 255, 0, 0 ) );
  cmap_grad2.addColor( DGtalBoard::Color( 255, 255, 10 ) );
  
  d=0;
  for( std::vector<Z3i::SCell>::iterator it=vectBdrySCell.begin(); 
       it!= vectBdrySCell.end(); it++){
     DGtalBoard::Color col= cmap_grad2(d);
     viewer.setFillColor(QColor(col.red(),col.green() ,col.blue(), 255));
     viewer<< *it;
    d++;
  }
  
  // On need once again to avoid superposition.
  viewer << DGtalQGLViewer::shiftSurfelVisu; 
  viewer.setFillColor(QColor(18, 200, 25, 255));
  viewer << aCell ;
  viewer << DGtalQGLViewer::updateDisplay;
    
  return application.exec();
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
