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
  this->m_StartTime = 0;
  this->m_StopTime = 0;
  this->m_Overlapped = false;
}

/** Constructor with default arguments */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::LocalAffineTransform(unsigned int parametersDimension):
  Superclass(parametersDimension)
{
  this->m_StartTime = 0;
  this->m_StopTime = 0;
  this->m_Overlapped = false;
}

/** Constructor with explicit arguments */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::LocalAffineTransform(const MatrixType & matrix,
                                                             const OutputVectorType & offset):
  Superclass(matrix, offset)
{
  this->m_StartTime = 0;
  this->m_StopTime = 0;
  this->m_Overlapped = false;
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
vnl_matrix<TScalarType>
LocalAffineTransform< TScalarType, NDimensions >::ComputeLogarithmMatrix(
  AffineTransformPointer affineTransform)
{
  vnl_matrix<TScalarType> homoMat(NDimensions+1, NDimensions+1);
  this->GetHomogeneousMatrix(homoMat, affineTransform);

  vnl_matrix<TScalarType> logMat = sdtools::GetLogarithm(homoMat);

  return logMat;
}

/** Compute the exponential of an affine transform */
template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformPointer
LocalAffineTransform< TScalarType, NDimensions >::ComputeExponentialTransform(
  vnl_matrix<TScalarType> velocity)
{
  AffineTransformPointer expAffine = AffineTransformType::New();

  vnl_matrix<TScalarType> expMat = sdtools::GetExponential(velocity);
  this->SetHomogeneousTransform(expAffine, expMat);

  return expAffine;
}

template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >::SetHomogeneousTransform(
  AffineTransformPointer &affineTransform,
  const vnl_matrix<TScalarType> &homoMatrix)
{
  MatrixType vmat;
  OutputVectorType voffset;
  TScalarType scale = homoMatrix[NDimensions][NDimensions];

  for (unsigned int i=0; i<NDimensions; i++)
    {
    for (unsigned int j=0; j<NDimensions; j++)
      {
        vmat[i][j] = homoMatrix[i][j] / scale;
      }
    voffset[i] = homoMatrix[i][NDimensions] / scale;
    }

  affineTransform->SetCenter(this->GetCenter());
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
vnl_matrix<TScalarType>
LocalAffineTransform< TScalarType, NDimensions >
::GetVelocityMatrix()
{
  if (this->m_VelocityMatrix.empty())
    {
    this->m_VelocityMatrix = this->ComputeLogarithmMatrix(this);
    }
  return this->m_VelocityMatrix;
}

template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >
::UpdatePartialVelocityMatrix()
{
  double timePeriod = (double)(m_StopTime - m_StartTime) / m_TimeStampMax;
  this->m_PartialVelocityMatrix = this->GetVelocityMatrix();
  this->m_PartialVelocityMatrix *= timePeriod;
}

template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::OutputVectorType
LocalAffineTransform< TScalarType, NDimensions >
::GetPartialVelocityAtPoint(const PointType &point)
{
  OutputVectorType velocity;
  
  for (unsigned int r=0; r<NDimensions; r++)
    {
    velocity[r] = this->m_PartialVelocityMatrix[r][NDimensions];
    for (unsigned int c=0; c<NDimensions; c++)
      {
      velocity[r] += this->m_PartialVelocityMatrix[r][c] * point[c];
      }
    }
  return velocity;
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
void
LocalAffineTransform< TScalarType, NDimensions >
::DilateFixedMaskImage(unsigned int radius)
  {
  typedef typename itk::BinaryBallStructuringElement<typename MaskImageType::PixelType,
    NDimensions> StructuringElementType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();
 
  typedef typename itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, 
    StructuringElementType> BinaryDilateImageFilterType;
 
  typename BinaryDilateImageFilterType::Pointer dilateFilter
          = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(this->m_FixedMaskImage);
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->SetDilateValue(1); //binary mask
  dilateFilter->Update();

  MaskImagePointer dilatedImage = dilateFilter->GetOutput();
  dilatedImage->DisconnectPipeline();
  this->m_FixedMaskImage = dilatedImage;
  }

template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformPointer 
LocalAffineTransform< TScalarType, NDimensions >
::GetPartialTransform(int startTime, int stopTime)
{
  AffineTransformPointer partialTransform = AffineTransformType::New();
  int duration = stopTime - startTime;
  if (duration == 0)
    {
    partialTransform->SetCenter(this->GetCenter());
    partialTransform->SetIdentity();
    }
  else if (duration == this->m_TimeStampMax)
    {
    partialTransform->SetCenter(this->GetCenter());
    partialTransform->SetMatrix(this->GetMatrix());
    partialTransform->SetOffset(this->GetOffset());
    }
  else if (duration == -this->m_TimeStampMax)
    {
    partialTransform->SetCenter(this->GetCenter());
    bool invertible = this->GetInverse(partialTransform);
    if (!invertible)
      {
      itkWarningMacro("This LocalAffineTransform is not invertible.");
      }
    }
  else
    {
    double factor = (double)duration / this->m_TimeStampMax;
    vnl_matrix<TScalarType> partialVelocity = this->GetVelocityMatrix();
    partialVelocity *= factor;
    partialTransform = this->ComputeExponentialTransform(partialVelocity);
    }

  return partialTransform;
}

/**
 * The existing method AffineTransform::Scale(factor) scales around the center.
 * Instead, we scale around (0,0) in ScaleMatrixOffset().
 */
template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformPointer 
LocalAffineTransform< TScalarType, NDimensions >
::ScaleMatrixOffset(AffineTransformPointer transform, double factor)
{
  AffineTransformPointer result = AffineTransformType::New();

  result->SetCenter(transform->GetCenter());

  typename AffineTransformType::MatrixType newMatrix = transform->GetMatrix();
  newMatrix *= factor;
  result->SetMatrix(newMatrix);
  
  typename AffineTransformType::OutputVectorType newOffset = transform->GetOffset();
  newOffset *= factor;
  result->SetOffset(newOffset);

  return result;
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
::WarpFixedMaskIntoPointSet(int timeStamp)
{
  AffineTransformPointer forwardTransform = this->GetPartialTransform(0, timeStamp);

  this->WarpMaskIntoPointSet(this->m_SamplePointSet,
    this->m_FixedMaskImage, forwardTransform);
}

template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::WarpMovingMaskIntoPointSet(int timeStamp)
{
  AffineTransformPointer backwardTransform = this->GetPartialTransform(this->m_TimeStampMax, timeStamp);

  this->WarpMaskIntoPointSet(this->m_SamplePointSet,
    this->m_MovingMaskImage, backwardTransform);
}

template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::ComputeSamplePointSet(int timeStamp)
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
