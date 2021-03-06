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
  this->m_TimeStampMax = 0;
}

/** Constructor with default arguments */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::LocalAffineTransform(unsigned int parametersDimension):
  Superclass(parametersDimension)
{
  this->m_StartTime = 0;
  this->m_StopTime = 0;
  this->m_Overlapped = false;
  this->m_TimeStampMax = 0;
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
  this->m_TimeStampMax = 0;
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

  this->m_TimerMatrixLogarithm.Start();
  vnl_matrix<TScalarType> logMat = sdtools::GetLogarithm(homoMat);
  this->m_TimerMatrixLogarithm.Stop();

  return logMat;
}

/** Compute the exponential of an affine transform */
template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformPointer
LocalAffineTransform< TScalarType, NDimensions >::ComputeExponentialTransform(
  vnl_matrix<TScalarType> velocity)
{
  AffineTransformPointer expAffine = AffineTransformType::New();

  this->m_TimerMatrixExponential.Start();
  vnl_matrix<TScalarType> expMat = sdtools::GetExponential(velocity);
  this->m_TimerMatrixExponential.Start();

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

/** Update the partial velocity matrix during current [m_StartTime, m_StopTime]. */
template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >
::UpdatePartialVelocityMatrix()
{
  if (this->m_TimeStampMax == 0)
    {
    itkExceptionMacro("LocalAffineTransform::m_TimeStampMax must be set before UpdatePartialVelocityMatrix is called.");
    return;
    }
  double timePeriod = (double)(m_StopTime - m_StartTime) / m_TimeStampMax;
  this->m_PartialVelocityMatrix = this->GetVelocityMatrix();
  this->m_PartialVelocityMatrix *= timePeriod;
}

/** Get the partial velocity vector at a point during current [m_StartTime, m_StopTime]. */ 
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
  this->m_TimerComputeFixedMaskImage.Start();

  typedef itk::SpatialObjectToImageFilter<TSpatialObject,MaskImageType> SpatialObjectToImageFilterType;
  typename SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
  imageFilter->SetInput(spatialObject);
  imageFilter->SetInsideValue(1);
  imageFilter->SetOutsideValue(0);
  imageFilter->SetSize(size);
  imageFilter->Update();
  this->m_FixedMaskImage = imageFilter->GetOutput();   

  this->m_TimerComputeFixedMaskImage.Stop();
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
  this->m_TimerComputeFixedMaskImage.Start();

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

  this->m_TimerComputeFixedMaskImage.Stop();
}

template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >
::DilateFixedMaskImage(unsigned int radius)
{
  m_TimerDilateFixedMaskImage.Start();

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

  m_TimerDilateFixedMaskImage.Stop();  
}

/** Compute the partial transform during [startTime, stopTime].
 *  It is computed by exp(t*log(T)) where t=(stopTime-startTime)/m_TimeStampMax,
 *  and T is this transform. */
template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::AffineTransformPointer 
LocalAffineTransform< TScalarType, NDimensions >
::GetPartialTransform(int startTime, int stopTime)
{
  if (this->m_TimeStampMax == 0)
    {
    itkExceptionMacro("LocalAffineTransform::m_TimeStampMax must be set before GetPartialTransform is called.");
    return NULL;
    }

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
      itkExceptionMacro("This LocalAffineTransform is not invertible.");
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

/** Warp a mask into a PointSet by a transform. */
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

/** Warp the fixed mask into a PointSet in a virtual domain at timeStamp. */
template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::WarpFixedMaskIntoPointSet(int timeStamp)
{
  AffineTransformPointer forwardTransform;
  if (timeStamp != 0)
    {
    forwardTransform = this->GetPartialTransform(0, timeStamp);
    }
  else
    {
    //it doesn't require m_TimeStampMax to be set at the beginning
    forwardTransform = AffineTransformType::New();
    forwardTransform->SetIdentity();
    }

  this->WarpMaskIntoPointSet(this->m_SamplePointSet,
    this->m_FixedMaskImage, forwardTransform);
}

/** Warp the moving mask into a PointSet in a virtual domain at timeStamp. */
template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::WarpMovingMaskIntoPointSet(int timeStamp)
{
  AffineTransformPointer backwardTransform;
  if (timeStamp != 0)
    {
    backwardTransform = this->GetPartialTransform(this->m_TimeStampMax, timeStamp);
    }
  else
    {
    //it doesn't require m_TimeStampMax to be set at the beginning
    backwardTransform = AffineTransformType::New();
    bool invertible = this->GetInverse(backwardTransform);
    if (!invertible)
      {
      itkExceptionMacro("This LocalAffineTransform is not invertible.");
      }
    }

  this->WarpMaskIntoPointSet(this->m_SamplePointSet,
    this->m_MovingMaskImage, backwardTransform);
}

template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::ComputeSamplePointSet(int timeStamp)
{
  m_TimerComputeSamplePointSet.Start();
  //fill m_SamplePointSet with points
  this->m_SamplePointSet = PointSetType::New();
  this->WarpFixedMaskIntoPointSet(timeStamp);
  this->WarpMovingMaskIntoPointSet(timeStamp);
  m_TimerComputeSamplePointSet.Stop();
}

/**
 * Transform the fixed mask into the moving mask. The moving mask
 * domain include both the fixed and moving domain.
 */
template< class TScalarType, unsigned int NDimensions >
void 
LocalAffineTransform< TScalarType, NDimensions >
::ComputeMovingMaskImage()
{
  this->m_TimerComputeMovingMaskImage.Start();

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
    this->m_TimerComputeMovingMaskImage.Stop();
    itkExceptionMacro("This LocalAffineTransform is not invertible.");
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

  this->m_TimerComputeMovingMaskImage.Stop();
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
