/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkLocalAffineTransform_hxx
#define __itkLocalAffineTransform_hxx

#include "itkNumericTraits.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl_sd_matrix_tools.h"
#include "itkLocalAffineTransform.h"

namespace itk
{
/** Constructor with default arguments */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::LocalAffineTransform():Superclass(ParametersDimension)
{}

/** Constructor with default arguments */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::LocalAffineTransform(unsigned int parametersDimension):
  Superclass(parametersDimension)
{}

/** Constructor with explicit arguments */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::LocalAffineTransform(const MatrixType & matrix,
                                                             const OutputVectorType & offset):
  Superclass(matrix, offset)
{}

/**  Destructor */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::
~LocalAffineTransform()
{
  return;
}

/** Compute principal logorithm of the homegeneous matrix */
template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >::ComputePrincipalLogorithm()
{
  this->m_VelocityAffineTransform = VelocityAffineTransformType::New();

  vnl_matrix<TScalarType> homoMat(NDimensions+1, NDimensions+1);
  MatrixType mat = this->GetMatrix();
  OutputVectorType offset = GetOffset();

  homoMat.fill(0.0);
  for (unsigned int i=0; i<NDimensions; i++)
    {
    for (unsigned int j=0; j<NDimensions; j++)
      {
      homoMat[i][j] = mat[i][j];
      }
    homoMat[i][NDimensions] = offset[i];
    }
  homoMat[NDimensions][NDimensions] = 1;

  std::cout << "homoMat = " << homoMat << std::endl;
  vnl_matrix<TScalarType> logMat = sdtools::GetLogarithm(homoMat);

  MatrixType vmat;
  OutputVectorType voffset;
  for (unsigned int i=0; i<NDimensions; i++)
    {
    for (unsigned int j=0; j<NDimensions; j++)
      {
        vmat[i][j] = logMat[i][j];
      }
    voffset[i] = logMat[i][NDimensions];
    }

  this->m_VelocityAffineTransform->SetMatrix(vmat);
  this->m_VelocityAffineTransform->SetOffset(voffset);

}

template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::VelocityAffineTransformType*
LocalAffineTransform< TScalarType, NDimensions >::GetVelocityAffineTransform()
{
  if (this->m_VelocityAffineTransform.IsNull())
    {
    this->ComputePrincipalLogorithm();
    }
  return this->m_VelocityAffineTransform.GetPointer();
}

template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >::AddFixedPoint(InputPointType point)
{
  if (this->m_FixedPointSet.IsNull())
    {
    this->m_FixedPointSet = PointSet::New();
    }
  PointSetType::PointIndentifier id = this->m_FixedPointSet->GetNumberOfPoints();
  this->m_FixedPointSet.SetPoint(id, point);
}

template< class TScalarType, unsigned int NDimensions >
template< class TSpatialObject > 
void 
LocalAffineTransform< TScalarType, NDimensions >
::ComputeFixedMaskImageFromSpatialObject(TSpatialObject *spatialObject, SizeType size)
{
  typedef itk::SpatialObjectToImageFilter<TSpatialObject,MaskImageType> SpatialObjectToImageFilterType;
  SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
  imageFilter->SetInput(spatialObject);
  imageFilter->SetInsideValue(1);
  imageFilter->SetOutsideValue(0);
  imageFilter->SetSize(size);
  imageFilter->Update();
  this->m_FixedMaskImage = imageFilter->GetOutput();   
}

template< class TScalarType, unsigned int NDimensions >
template< class TSpatialObject > 
void 
LocalAffineTransform< TScalarType, NDimensions >
::ComputeFixedMaskImageFromSpatialObject(TSpatialObject *spatialObject,
                                PointType origin,
                                DirectionType direction,
                                SpacingType spacing,
                                SizeType size)
{
  typedef itk::SpatialObjectToImageFilter<TSpatialObject,MaskImageType> SpatialObjectToImageFilterType;
  SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
  imageFilter->SetInput(spatialObject);
  imageFilter->SetInsideValue(1);
  imageFilter->SetOutsideValue(0);

  imageFilter->SetOrigin(origin);
  imageFilter->SetSpacing(spacing);
  imageFilter->SetDirection(direction);
  imageFilter->SetSize(size);

  imageFilter->Update();
  this->m_FixedMaskImage = imageFilter->GetOutput();
}

template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::AddMaskToPointSet(PointSetPointer &pointSet, const MaskImagePointer &mask)
{
  typedef ImageRegionIteratorWithIndex< MaskImageType > IteratorType;
  IteratorType it( mask, mask->GetLargestPossibleRegion() );

  //fill the moving mask with values
  PointSetType::PointIdentifier pointId = pointSet->GetNumberOfPoints();
  PointType point;

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
      if (it.Get() != 0) //foreground
        {
        mask->TransformIndexToPhysicalPoint(it.GetIndex(), point);
        pointSet->SetPoint(pointId++, point);
      }
    }
}

/**
 * Warp the fixed mask to the moving mask. The moving mask
 * uses the meta information from the fixed mask except its
 * region. Instead, the region will be expanded to contain
 * the warped mask.
 */
template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::ComputeMovingMaskImageAndDenseFixedPointSet()
{
  ContinuousIndexType mappedCorner, minIndex, maxIndex;

  //compute minIndex and maxIndex of the warped mask
  RegionType fixedRegion = this->m_FixedMaskImage->GetLargestPossibleRegion();
  minIndex = fixedRegion.GetIndex();
  maxIndex = fixedRegion.GetUpperIndex();

  typedef typename itk::VectorContainer<int, IndexType> IndexContainerType;
  typename IndexContainerType::Pointer corners = 
    PicslImageHelper::GetCorners<RegionType>(fixedRegion);
  IndexContainerType::ConstIterator corner;
  for ( corner = corners->Begin(); corner != corners->End(); corner++)
    {
    IndexType cornerIndex = corner.Value();
    PointType cornerPoint, mappedPoint;
    //map the corner index into the physical point
    this->m_FixedMaskImage->TransformIndexToPhysicalPoint(cornerIndex, cornerPoint);
    //transform the point by this local affine transform
    mappedPoint = this->TransformPoint(cornerPoint);
    //map the transformed point to index
    this->m_FixedMaskImage->TransformPhysicalPointToContinuousIndex(mappedPoint, mappedCorner);

    CopyWithMin<ContinuousIndexType>(minIndex, mappedCorner);
    CopyWithMax<ContinuousIndexType>(maxIndex, mappedCorner);
    }
  
  //allocate the moving mask with the possible region
  IndexType index1, index2; 
  index1.CopyWithRound(minIndex);
  index2.CopyWithRound(maxIndex);
  RegionType movingRegion;
  movingRegion.SetIndex(index1);
  movingRegion.SetUpperIndex(index2);

  MaskImagePointer movingMask = MaskImageType::New();
  movingMask->CopyInformation(this->m_FixedMaskImage);
  movingMask->SetRegions(movingRegion);
  movingMask->Allocate();
  
  //compute the inverse transform
  Superclass::Pointer inverse = Superclass::New();
  bool insideImage, invertible = this->GetInverse(inverse);
  if (!invertible)
    {
    itkWarningMacro("This LocalAffineTransform is not invertible.");
    return;
    }

  PointType point1, point2;
  MaskImageType::PixelType pixel;

  IndexType minMovingMaskIndex, maxMovingMaskIndex;
  minMovingMaskIndex = fixedRegion.GetIndex();
  maxMovingMaskIndex = fixedRegion.GetUpperIndex();

  //fill m_DenseFixedPointSet with points
  this->m_DenseFixedPointSet = PointSetType::New();
  this->AddMaskToPointSet(this->m_DenseFixedPointSet, this->m_FixedMaskImage);
  int pointId = this->m_DenseFixedPointSet->GetNumberOfPoints();

  //fill the moving mask with values
  typedef ImageRegionIteratorWithIndex< MaskImageType > IteratorType;
  IteratorType it( movingMask, movingMask->GetLargestPossibleRegion() );
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index1 = it.GetIndex();
    movingMask->TransformIndexToPhysicalPoint(index1, point1);
    point2 = inverse->TransformPoint(point1);

    insideImage = this->m_FixedMaskImage->TransformPhysicalPointToIndex(point2, index2);
    pixel = 0;
    if (insideImage)
      {
      pixel = this->m_FixedMaskImage->GetPixel(index2);
      if (pixel != 0) //foreground
        {
        //update min and max indices
        CopyWithMin<IndexType>(minMovingMaskIndex, index1);
        CopyWithMax<IndexType>(maxMovingMaskIndex, index1);
        this->m_DenseFixedPointSet->SetPoint(pointId++, point2);
        }
      }
    it.Set(pixel);
    } //for iterator of movingMask

  RegionType extractRegion;
  extractRegion.SetIndex(minMovingMaskIndex);
  extractRegion.SetUpperIndex(maxMovingMaskIndex);

  typedef itk::ExtractImageFilter< MaskImageType, MaskImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetExtractionRegion(extractRegion);
  filter->SetInput(movingMask);
#if ITK_VERSION_MAJOR >= 4
  filter->SetDirectionCollapseToIdentity(); // This is required.
#endif
  filter->Update();

  this->m_MovingMaskImage = filter->GetOutput();
}

template< class TScalarType, unsigned int NDimensions >
template< class TVector>
static void
LocalAffineTransform< TScalarType, NDimensions >
::CopyWithMax(TVector &maxVec, const TVector &newVec)
{
  for (unsigned int d=0; d<NDimensions; d++)
    {
    if (maxVec[d] < newVec[d])
      {
      maxVec[d] = newVec[d];
      }
    }
}

template< class TScalarType, unsigned int NDimensions >
template< class TVector>
void
LocalAffineTransform< TScalarType, NDimensions >
::CopyWithMin(TVector &minVec, const TVector &newVec)
{
  for (unsigned int d=0; d<NDimensions; d++)
    {
    if (minVec[d] > newVec[d])
      {
      minVec[d] = newVec[d];
      }
    }
}

/** Print self */
template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // namespace

#endif
