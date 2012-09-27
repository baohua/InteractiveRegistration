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
  this->m_FixedPointSet.push_back(point);
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

/**
 * Warp the fixed mask to the moving mask. The moving mask
 * uses the meta information from the fixed mask except its
 * region. Instead, the region will be expanded to contain
 * the warped mask.
 */
template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::ComputeMovingMaskImageFromFixedMaskImage()
{
  ContinuousIndexType mappedCorner, minIndex, maxIndex;
  PointType point, mappedPoint;

  //compute minIndex and maxIndex of the warped mask
  RegionType fixedRegion = this->m_FixedMaskImage->GetLargestPossibleRegion();
  minIndex = fixedRegion.GetIndex();
  maxIndex = fixedRegion.GetUpperIndex();

  int cornerNum = 1 << NDimensions;
  IndexType corner, firstCorner = fixedRegion.GetIndex();
  SizeType size = fixedRegion.GetSize();
  for(int i=0; i<cornerNum; i++)
    {
    int bit;
    for (int d=0; d<NDimensions; d++)
      {
      bit = (int) (( i & (1 << d) ) != 0); // 0 or 1
      corner[d] = firstCorner[d] + bit * (size[d] - 1);
      }
    //map the corner index into the physical point
    this->m_FixedMaskImage->TransformIndexToPhysicalPoint(corner, point);
    //transform the point by this local affine transform
    mappedPoint = this->TransformPoint(point);
    //map the transformed point to index
    this->m_FixedMaskImage->TransformPhysicalPointToContinuousIndex(mappedPoint, mappedCorner);

    CopyWithMin<ContinuousIndexType>(minIndex, mappedCorner);
    CopyWithMax<ContinuousIndexType>(maxIndex, mappedCorner);
    } //for each corner 
  
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

  typedef ImageRegionIteratorWithIndex< MaskImageType > IteratorType;
  IteratorType it( movingMask, movingMask->GetLargestPossibleRegion() );

  //fill the moving mask with values
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index1 = it.GetIndex();
    movingMask->TransformIndexToPhysicalPoint(index1, point1);
    point2 = inverse->TransformPoint(point1);

    insideImage = this->m_FixedMaskImage->TransformPhysicalPointToIndex(point2, index2);
    if (insideImage)
      {
      pixel = this->m_FixedMaskImage->GetPixel(index2);
      it.Set(pixel);

      if (pixel > 0) //foreground
        {
        //update min and max indices
        CopyWithMin<IndexType>(minMovingMaskIndex, index1);
        CopyWithMax<IndexType>(maxMovingMaskIndex, index1);
        }
      }
    else
      {
      it.Set(0);
      }
    } //for iterator of movingMask

  RegionType trimmedRegion;
  trimmedRegion.SetIndex(minMovingMaskIndex);
  trimmedRegion.SetUpperIndex(maxMovingMaskIndex);

  typedef itk::ExtractImageFilter< MaskImageType, MaskImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetExtractionRegion(trimmedRegion);
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
