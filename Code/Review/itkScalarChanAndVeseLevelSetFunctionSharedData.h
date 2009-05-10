/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkScalarChanAndVeseLevelSetFunctionSharedData.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkScalarChanAndVeseLevelSetFunctionSharedData_h
#define __itkScalarChanAndVeseLevelSetFunctionSharedData_h

#include "itkLightObject.h"

#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkEuclideanDistance.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{

/** \class ScalarChanAndVeseLevelSetFunctionSharedData
 *
 * \brief Helper class used to share data in the ScalarChanAndVeseLevelSetFunction.
 *
 * This class holds cache data used during the computation of the LevelSet updates.
 *
 * Based on the paper:
 *
 *        "An active contour model without edges"
 *         T. Chan and L. Vese. 
 *         In Scale-Space Theories in Computer Vision, pages 141–151, 1999.
 * 
 * \author Mosaliganti K., Smith B., Gelas A., Gouaillard A., Megason S.
 *
 *  This code was taken from the Insight Journal paper:
 *
 *      "Cell Tracking using Coupled Active Surfaces for Nuclei and Membranes"
 *      http://www.insight-journal.org/browse/publication/642
 *      http://hdl.handle.net/10380/3055
 *
 *  That is based on the papers:
 *
 *      "Level Set Segmentation: Active Contours without edge"
 *      http://www.insight-journal.org/browse/publication/322
 *      http://hdl.handle.net/1926/1532
 *
 *      and
 *
 *      "Level set segmentation using coupled active surfaces"
 *      http://www.insight-journal.org/browse/publication/323
 *      http://hdl.handle.net/1926/1533
 *
 *
 */
template < class TInputImage, class TFeatureImage >
class ScalarChanAndVeseLevelSetFunctionSharedData : public LightObject
{
public:

  typedef ScalarChanAndVeseLevelSetFunctionSharedData       Self;
  typedef LightObject                                       Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  itkStaticConstMacro( ImageDimension, unsigned int, TFeatureImage::ImageDimension );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  itkTypeMacro(ScalarChanAndVeseLevelSetFunctionSharedData, LightObject);

  typedef TInputImage                                   InputImageType;
  typedef typename InputImageType::Pointer              InputImagePointer;
  typedef typename InputImageType::ConstPointer         InputImageConstPointer;
  typedef typename InputImageType::PixelType            InputPixelType;
  typedef typename InputImageType::RegionType           InputRegionType;
  typedef typename InputImageType::SizeType             InputSizeType;
  typedef typename InputSizeType::SizeValueType         InputSizeValueType;
  typedef typename InputImageType::SpacingType          InputSpacingType;
  typedef typename InputImageType::IndexType            InputIndexType;
  typedef typename InputIndexType::IndexValueType       InputIndexValueType;
  typedef typename InputImageType::PointType            InputPointType;

  typedef TFeatureImage                                 FeatureImageType;
  typedef typename FeatureImageType::Pointer            FeatureImagePointer;
  typedef typename FeatureImageType::ConstPointer       FeatureImageConstPointer;
  typedef typename FeatureImageType::PixelType          FeaturePixelType;
  typedef typename FeatureImageType::RegionType         FeatureRegionType;
  typedef typename FeatureImageType::SizeType           FeatureSizeType;
  typedef typename FeatureSizeType::SizeValueType       FeatureSizeValueType;
  typedef typename FeatureImageType::SpacingType        FeatureSpacingType;
  typedef typename FeatureImageType::IndexType          FeatureIndexType;
  typedef typename FeatureImageType::PointType          FeaturePointType;

  typedef std::list< unsigned int >                     ListPixelType;
  typedef Image< ListPixelType, ImageDimension >        ListImageType;
  typedef typename ListImageType::Pointer               ListImagePointer;
  typedef typename ListImageType::ConstPointer          ListImageConstPointer;
  typedef typename ListImageType::RegionType            ListRegionType;
  typedef typename ListImageType::SizeType              ListSizeType;
  typedef typename ListSizeType::SizeValueType          ListSizeValueType;
  typedef typename ListImageType::SpacingType           ListSpacingType;
  typedef typename ListImageType::IndexType             ListIndexType;
  typedef typename ListIndexType::IndexValueType        ListIndexValueType;
  typedef typename ListImageType::PointType             ListPointType;
  typedef ImageRegionIteratorWithIndex< ListImageType > ListIteratorType;

  typedef Vector< float, ImageDimension >                     CentroidVectorType;
  typedef itk::Statistics::ListSample< CentroidVectorType >   SampleType;
  typedef itk::Statistics::KdTreeGenerator< SampleType >      TreeGeneratorType;
  typedef typename TreeGeneratorType::Pointer                 TreePointer;
  typedef typename TreeGeneratorType::KdTreeType              TreeType;
  typedef typename TreeType::Pointer                          KdTreePointer;

  void SetFunctionCount( const unsigned int& n )
    {
    this->m_FunctionCount = n;

    this->m_CVals.resize( n, 0.0 );
    this->m_CDens.resize( n, 0.0 );
    this->m_CNums.resize( n, 0.0 );

    this->m_CB.resize( n, 0.0 );
    this->m_CBDen.resize( n, 0.0 );
    this->m_CBNum.resize( n, 0.0 );

    this->m_HVals.resize( n, 0 );
    this->m_Start.resize( n );
    this->m_End.resize( n );
    }

  void CreateHVals( const unsigned int& j,
    const InputSpacingType& spacing,
    const InputPointType& origin,
    const InputRegionType& region )
    {
    this->m_HVals[j] = InputImageType::New();
    this->m_HVals[j]->SetRegions( region );
    this->m_HVals[j]->Allocate();
    this->m_HVals[j]->SetOrigin( origin );
    this->m_HVals[j]->SetSpacing( spacing );
    this->m_HVals[j]->FillBuffer( 0 );

    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      this->m_Start[j][i] = static_cast< InputIndexValueType >( origin[i]/spacing[i] );
      this->m_End[j][i] = this->m_Start[j][i] + static_cast< InputIndexValueType >( region.GetSize()[i] ) - 1;
      }
    }

  void SetKdTree( KdTreePointer kdtree )
    {
    this->m_KdTree = kdtree;
    }

  template< class TIndex >
  bool VerifyInsideRegion( const unsigned int& i, const TIndex& featureIndex )
    {
    typedef typename TIndex::IndexValueType TIndexValueType;
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      if(  (featureIndex[j] < static_cast< TIndexValueType >(this->m_Start[i][j]) )
        || (featureIndex[j] > static_cast< TIndexValueType >(this->m_End[i][j]))  )
        {
        return false;
        }
      }
    return true;
    }

  InputIndexType GetIndex( const unsigned int& j, const FeatureIndexType& featureIndex )
    {
    InputIndexType index;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      index[i] = featureIndex[i] - static_cast< InputIndexValueType >( this->m_Start[j][i] );
      }

    return index;
    }

  FeatureIndexType GetFeatureIndex( const unsigned int& j,
    const InputIndexType& inputIndex )
    {
    FeatureIndexType index;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      index[i] = inputIndex[i] +
        static_cast< InputIndexValueType >( this->m_Start[j][i] );

    return index;
    }


  void AllocateListImage( FeatureRegionType region, FeatureSpacingType spacing )
    {
    // FIXME: Are we missing Origin and Orientation ?
    this->m_LImage = ListImageType::New();
    this->m_LImage->SetRegions( region );
    this->m_LImage->Allocate();
    this->m_LImage->SetSpacing( spacing );
    }

  void PopulateListImage()
    {
    ListSpacingType spacing = this->m_LImage->GetSpacing();
    ListIteratorType lIt( this->m_LImage, this->m_LImage->GetLargestPossibleRegion() );

    if ( m_KdTree )
      {
      for(lIt.GoToBegin(); !lIt.IsAtEnd(); ++lIt )
        {
        ListIndexType ind = lIt.GetIndex();

        unsigned int numberOfNeighbors = 6;
        float queryPoint[ImageDimension];
        for( unsigned int i = 0; i < ImageDimension; i++ )
          queryPoint[i] = ind[i]*spacing[i];

        typename TreeType::InstanceIdentifierVectorType neighbors;
        this->m_KdTree->Search( queryPoint, numberOfNeighbors, neighbors );

        ListPixelType L;
        for( unsigned int i = 0; i < numberOfNeighbors; i++ )
          {
          if( VerifyInsideRegion( neighbors[i], ind ) )
            {
            L.push_back( neighbors[i] );
            }
          }
        lIt.Set( L );
        }
      }
    else
      {
      for(lIt.GoToBegin(); !lIt.IsAtEnd(); ++lIt )
        {
        ListIndexType ind = lIt.GetIndex();
        ListPixelType L;
        for( unsigned int i = 0; i < this->m_FunctionCount; i++ )
          {
          if( VerifyInsideRegion( i, ind ) )
            {
            L.push_back( i );
            }
          }
        lIt.Set( L );
        }
      }
    }

  std::vector< double >             m_CB;
  std::vector< double >             m_CVals;
  std::vector< double >             m_CNums;
  std::vector< double >             m_CDens;
  std::vector< double >             m_CBNum;
  std::vector< double >             m_CBDen;

  unsigned int                      m_FunctionCount;
  std::vector< InputImagePointer >  m_HVals;
  std::vector< InputIndexType >     m_Start;
  std::vector< InputIndexType >     m_End;
  ListImagePointer                  m_LImage;
  KdTreePointer                     m_KdTree;

protected:
  ScalarChanAndVeseLevelSetFunctionSharedData() {}
  ~ScalarChanAndVeseLevelSetFunctionSharedData(){}

private:
  ScalarChanAndVeseLevelSetFunctionSharedData(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} //end namespace itk

#endif
