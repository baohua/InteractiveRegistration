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
{
  this->m_StartTime = 0.0;
  this->m_StopTime = 0.0;
  this->m_TimePeriod = 0.0;
}

/** Constructor with default arguments */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::LocalAffineTransform(unsigned int parametersDimension):
  Superclass(parametersDimension)
{
  this->m_StartTime = 0.0;
  this->m_StopTime = 0.0;
  this->m_TimePeriod = 0.0;
}

/** Constructor with explicit arguments */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::LocalAffineTransform(const MatrixType & matrix,
                                                             const OutputVectorType & offset):
  Superclass(matrix, offset)
{
  this->m_StartTime = 0.0;
  this->m_StopTime = 0.0;
  this->m_TimePeriod = 0.0;
}

/**  Destructor */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::
~LocalAffineTransform()
{
  return;
}

/** Compute the principal logarithm of an affine transform */
template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformPointer
LocalAffineTransform< TScalarType, NDimensions >::ComputeLogarithmTransform(
  AffineTransformPointer affineTransform)
{
  AffineTransformPointer logAffine = AffineTransformType::New();

  vnl_matrix<TScalarType> homoMat(NDimensions+1, NDimensions+1);
  this->GetHomogeneousMatrix(homoMat, affineTransform);

  vnl_matrix<TScalarType> logMat = sdtools::GetLogarithm(homoMat);
  this->SetHomogeneousTransform(logAffine, logMat, affineTransform->GetCenter());

  return logAffine;
}

/** Compute the exponential of an affine transform */
template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformPointer
LocalAffineTransform< TScalarType, NDimensions >::ComputeExponentialTransform(
  AffineTransformPointer affineTransform)
{
  AffineTransformPointer expAffine = AffineTransformType::New();

  vnl_matrix<TScalarType> homoMat(NDimensions+1, NDimensions+1);
  this->GetHomogeneousMatrix(homoMat, affineTransform);

  vnl_matrix<TScalarType> expMat = sdtools::GetExponential(homoMat);
  this->SetHomogeneousTransform(expAffine, expMat, affineTransform->GetCenter());

  return expAffine;
}

template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >::SetHomogeneousTransform(
  AffineTransformPointer &affineTransform,
  const vnl_matrix<TScalarType> &homoMatrix, InputPointType center)
{
  MatrixType vmat;
  OutputVectorType voffset;
  for (unsigned int i=0; i<NDimensions; i++)
    {
    for (unsigned int j=0; j<NDimensions; j++)
      {
        vmat[i][j] = homoMatrix[i][j];
      }
    voffset[i] = homoMatrix[i][NDimensions];
    }

  affineTransform->SetCenter(center);
  affineTransform->SetMatrix(vmat);
  affineTransform->SetOffset(voffset);
}

template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >::GetHomogeneousMatrix(
  vnl_matrix<TScalarType> &homoMatrix, 
  const AffineTransformPointer &affineTransform)
{
  MatrixType mat = affineTransform->GetMatrix();
  OutputVectorType offset = affineTransform->GetOffset();

  homoMatrix.fill(0.0);
  for (unsigned int i=0; i<NDimensions; i++)
    {
    for (unsigned int j=0; j<NDimensions; j++)
      {
      homoMatrix[i][j] = mat[i][j];
      }
    homoMatrix[i][NDimensions] = offset[i];
    }
  homoMatrix[NDimensions][NDimensions] = 1;  
}

template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformType*
LocalAffineTransform< TScalarType, NDimensions >
::GetVelocityAffineTransform()
{
  if (this->m_VelocityAffineTransform.IsNull())
    {
    this->m_VelocityAffineTransform = this->ComputeLogarithmTransform(this);
    }
  return this->m_VelocityAffineTransform.GetPointer();
}

template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformType*
LocalAffineTransform< TScalarType, NDimensions >
::GetPartialVelocityAffineTransform(double timePeriod)
{
  if (this->m_PartialVelocityAffineTransform.IsNull()
      || this->m_TimePeriod != timePeriod)
    {
    this->m_TimePeriod = timePeriod;
    this->m_PartialVelocityAffineTransform = this->ScaleMatrixOffset(
      this->GetVelocityAffineTransform(), timePeriod);
    }

  return this->m_PartialVelocityAffineTransform.GetPointer();
}

template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >::AddFixedPoint(InputPointType point)
{
  if (this->m_FixedPointSet.IsNull())
    {
    this->m_FixedPointSet = PointSetType::New();
    }
  typename PointSetType::PointIndentifier id = this->m_FixedPointSet->GetNumberOfPoints();
  this->m_FixedPointSet.SetPoint(id, point);
}

template< class TScalarType, unsigned int NDimensions >
template< class TSpatialObject > 
void 
LocalAffineTransform< TScalarType, NDimensions >
::ComputeFixedMaskImageFromSpatialObject(TSpatialObject *spatialObject, SizeType size)
{
  typedef itk::SpatialObjectToImageFilter<TSpatialObject,MaskImageType> SpatialObjectToImageFilterType;
  typename SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
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
  typename SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
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
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformPointer 
LocalAffineTransform< TScalarType, NDimensions >
::GetPartialTransform(double factor)
{
  AffineTransformPointer partialTransform = AffineTransform::New();
  if (factor == 0)
    {
    partialTransform->SetCenter(this->GetCenter());
    partialTransform->SetIdentity();
    return partialTransform;
    }
  else if (factor == 1)
    {
    partialTransform->SetCenter(this->GetCenter());
    partialTransform->SetMatrix(this->GetMatrix());
    partialTransform->SetOffset(this->GetOffset());
    return partialTransform;
    }
  else if (factor == -1)
    {
    partialTransform->SetCenter(this->GetCenter());
    bool invertible = this->GetInverse(partialTransform);
    if (!invertible)
      {
      itkWarningMacro("This LocalAffineTransform is not invertible.");
      }
    return partialTransform;
    }

  AffineTransformPointer partialVelocity = this->ScaleMatrixOffset(
    this->GetVelocityAffineTransform(), factor);

  partialTransform = this->ComputeExponentialTransform(partialVelocity);

  return partialTransform;
}

/**
 * The existing method AffineTransform::Scale(factor) scales around the center.
 * Instead, we scale around (0,0) in ScaleMatrixOffset().
 */
template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformPointer 
LocalAffineTransform< TScalarType, NDimensions >
::ScaleMatrixOffset(AffineTransformPointer velocity, double factor)
{
  AffineTransformPointer partialVelocity = AffineTransformType::New();

  partialVelocity->SetCenter(velocity->GetCenter());

  AffineTransformType::MatrixType newMatrix = velocity->GetMatrix();
  newMatrix *= factor;
  partialVelocity->SetMatrix(newMatrix);
  
  AffineTransformType::OutputVectorType newOffset = velocity->GetOffset();
  newOffset *= factor;
  partialVelocity->SetOffset(newOffset);

  return partialVelocity;
}

template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformPointer 
LocalAffineTransform< TScalarType, NDimensions >
::GetResidualTransform(double portion)
{
  AffineTransformPointer partialVelocity = this->GetPartialTransform(1-portion);

  return partialVelocity;
}

template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::WarpMaskIntoPointSet(PointSetPointer &pointSet, const MaskImagePointer &mask,
                       const AffineTransformPointer &transform)
{
  typedef ImageRegionIteratorWithIndex< MaskImageType > IteratorType;
  IteratorType it( mask, mask->GetLargestPossibleRegion() );

  //fill the moving mask with values
  typename PointSetType::PointIdentifier pointId = pointSet->GetNumberOfPoints();
  PointType point, point2;

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
      if (it.Get() != 0) //foreground
        {
        mask->TransformIndexToPhysicalPoint(it.GetIndex(), point);
        point2 = transform->TransformPoint(point); 
        pointSet->SetPoint(pointId++, point2);
      }
    }
}

template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::WarpFixedMaskIntoPointSet(double timeStamp)
{
  AffineTransformPointer forwardTransform = this->GetPartialTransform(timeStamp);

  this->WarpMaskIntoPointSet(this->m_SamplePointSet,
    this->m_FixedMaskImage, forwardTransform);
}

template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::WarpMovingMaskIntoPointSet(double timeStamp)
{
  AffineTransformPointer backwardTransform = this->GetPartialTransform(timeStamp-1);

  this->WarpMaskIntoPointSet(this->m_SamplePointSet,
    this->m_MovingMaskImage, backwardTransform);
}

template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::ComputeSamplePointSet(double timeStamp)
{
  //fill m_SamplePointSet with points
  this->m_SamplePointSet = PointSetType::New();
  this->WarpFixedMaskIntoPointSet(timeStamp);
  this->WarpMovingMaskIntoPointSet(timeStamp);
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
::ComputeMovingMaskImage()
{
  ContinuousIndexType mappedCorner, minIndex, maxIndex;

  //compute minIndex and maxIndex of the warped mask
  RegionType fixedRegion = this->m_FixedMaskImage->GetLargestPossibleRegion();
  minIndex = fixedRegion.GetIndex();
  maxIndex = fixedRegion.GetUpperIndex();

  typedef typename itk::VectorContainer<int, IndexType> IndexContainerType;
  typename IndexContainerType::Pointer corners = 
    PicslImageHelper::GetCorners<RegionType>(fixedRegion);
  typename IndexContainerType::ConstIterator corner;
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

    PicslImageHelper::CopyWithMin<ContinuousIndexType>(minIndex, mappedCorner);
    PicslImageHelper::CopyWithMax<ContinuousIndexType>(maxIndex, mappedCorner);
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
  typename Superclass::Pointer inverse = Superclass::New();
  bool insideImage, invertible = this->GetInverse(inverse);
  if (!invertible)
    {
    itkWarningMacro("This LocalAffineTransform is not invertible.");
    return;
    }

  PointType point1, point2;
  typename MaskImageType::PixelType pixel;

  IndexType minMovingMaskIndex, maxMovingMaskIndex;
  minMovingMaskIndex = fixedRegion.GetIndex();
  maxMovingMaskIndex = fixedRegion.GetUpperIndex();

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
        PicslImageHelper::CopyWithMin<IndexType>(minMovingMaskIndex, index1);
        PicslImageHelper::CopyWithMax<IndexType>(maxMovingMaskIndex, index1);
        }
      }
    it.Set(pixel);
    } //for iterator of movingMask

  RegionType extractRegion;
  extractRegion.SetIndex(minMovingMaskIndex);
  extractRegion.SetUpperIndex(maxMovingMaskIndex);

  typedef itk::ExtractImageFilter< MaskImageType, MaskImageType > FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetExtractionRegion(extractRegion);
  filter->SetInput(movingMask);
#if ITK_VERSION_MAJOR >= 4
  filter->SetDirectionCollapseToIdentity(); // This is required.
#endif
  filter->Update();

  this->m_MovingMaskImage = filter->GetOutput();
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
