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
#ifndef __itkPolyAffineTransform_hxx
#define __itkPolyAffineTransform_hxx

#include "itkNumericTraits.h"
#include "itkPolyAffineTransform.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "itkMath.h"
#include "itkCrossHelper.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
// Constructor with default arguments
template< class TScalarType,
          unsigned int NDimensions >
PolyAffineTransform< TScalarType, NDimensions >
::PolyAffineTransform():
  Superclass(0)
{
  this->SetTimeStampLog(8); //by default 2^8 time stamps
  this->m_TimeStampAutomaticFlag = true; //true to overwrite SetTimeStampLog()

  this->m_PadBoundary = 40;
  this->m_PadTrajectory = 2;

  this->m_DecayRateOfBoundary = 1*1;
  this->m_DecayRateOfTrajectory = 3*3;

  this->m_StopAllTrajectoriesAtOverlap = false;
}

// Destructor
template< class TScalarType, unsigned int NDimensions >
PolyAffineTransform< TScalarType, NDimensions >
::~PolyAffineTransform()
{
  return;
}

// Print self
template< class TScalarType, unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  unsigned int t;

  os << indent << "PolyAffineTransform: contains "
     << m_LocalAffineTransformVector.size() << " transforms " << std::endl;

  for ( t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    os << indent << "Transforms[" << t << "] = ";
    os << indent.GetNextIndent();
    //m_LocalAffineTransformVector[t]->PrintSelf(os, indent.GetNextIndent());
    os << m_LocalAffineTransformVector[t]->GetTransformTypeAsString();
    os << std::endl;
    os << indent.GetNextIndent() << "NumberOfParameters ";
    os << m_LocalAffineTransformVector[t]->GetNumberOfParameters();
    os << std::endl;
    }
}
// Constructor with explicit arguments
template< class TScalarType, unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::SetIdentity(void)
{
  for ( unsigned int t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    m_LocalAffineTransformVector[t]->SetIdentity();
    }
  this->GetParameters();
  this->GetFixedParameters();
  this->Modified();
}

// Transform a point
template< class TScalarType, unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::OutputPointType
PolyAffineTransform< TScalarType, NDimensions >
::TransformPoint(const InputPointType & point) const
{
  if (m_LocalAffineTransformVector.size() <= 0)
    {
    return point;
    }

  unsigned int i;
  OutputPointType outPoint, sumPoint;
  double squaredDistance, weight, sumWeight;

  sumPoint.Fill(0.0);
  sumWeight = 0.0;

  for ( i = 0; i < m_LocalAffineTransformVector.size(); i++ )
    {
    outPoint = m_LocalAffineTransformVector[i]->TransformPoint( point );

    squaredDistance = point.SquaredEuclideanDistanceTo(m_LocalAffineTransformVector[i]->GetCenter());
    weight = vcl_exp( - squaredDistance );

    sumWeight = sumWeight + weight;
    for ( unsigned d = 0; d < NDimensions; d++ )
      {
      sumPoint[d] += outPoint[d] * weight;
      }
    }

  for ( unsigned d = 0; d < NDimensions; d++ )
    {
    sumPoint[d] /= sumWeight;
    }

  return sumPoint;
}

// Transform a vector
template< class TScalarType, unsigned int NDimensions >
typename PolyAffineTransform< TScalarType,
                                    NDimensions >::OutputVectorType
PolyAffineTransform< TScalarType, NDimensions >
::TransformVector(const InputVectorType & vect) const
{
  itkExceptionMacro("TransformVector not yet implemented.");
  OutputVectorType output;
  return output;
}

// Transform a vnl_vector_fixed
template< class TScalarType, unsigned int NDimensions >
typename PolyAffineTransform< TScalarType,
                                    NDimensions >::OutputVnlVectorType
PolyAffineTransform< TScalarType, NDimensions >
::TransformVector(const InputVnlVectorType & vect) const
{
  itkExceptionMacro("TransformVector not yet implemented.");
  OutputVnlVectorType output;
  return output;
}

// Transform a CovariantVector
template< class TScalarType, unsigned int NDimensions >
typename PolyAffineTransform< TScalarType,
                                    NDimensions >::OutputCovariantVectorType
PolyAffineTransform< TScalarType, NDimensions >
::TransformCovariantVector(const InputCovariantVectorType & vec) const
{
  itkExceptionMacro("TransformCovariantVector not yet implemented.");
  OutputCovariantVectorType result;     // Converted vector
  return result;
}

// return an inverse transformation
template< class TScalarType,
          unsigned int NDimensions >
bool
PolyAffineTransform< TScalarType, NDimensions >
::GetInverse(Self *inverse) const
{
  if ( !inverse )
    {
    return false;
    }

  unsigned int i;

  for ( i = 0; i < m_LocalAffineTransformVector.size(); i++ )
    {
      if ( !m_LocalAffineTransformVector[i]->GetInverse( inverse->GetLocalAffineTransformVector()[i].GetPointer() ) )
        {
        return false;
        }
    }

  return true;
}

// Return an inverse of this transform
template< class TScalarType, unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::InverseTransformBasePointer
PolyAffineTransform< TScalarType, NDimensions >
::GetInverseTransform() const
{
  Pointer inv = New();

  return GetInverse(inv) ? inv.GetPointer() : NULL;
}

// Get fixed parameters
template< class TScalarType, unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::SetFixedParameters(const ParametersType & fp)
{
  unsigned int par = 0;

  //Save parameters. Needed for proper operation of TransformUpdateParameters.
  if( &fp != &(this->m_FixedParameters) )
    {
    this->m_FixedParameters = fp;
    }

  unsigned int paramSize;

  for ( unsigned int t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    //assuming each transform object has m_FixedParameters set initially
    LocalAffineParametersType param = m_LocalAffineTransformVector[t]->GetFixedParameters();
    paramSize = param.GetSize();

    for ( unsigned int par1 = 0; par1 < paramSize; par1++ )
      {
      param[par1] = this->m_FixedParameters[par];
      par++;
      }

    m_LocalAffineTransformVector[t]->SetFixedParameters(param);
    }

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();
}

/** Get the Fixed Parameters. */
template< class TScalarType, unsigned int NDimensions >
const typename PolyAffineTransform< TScalarType,
                                          NDimensions >::ParametersType &
PolyAffineTransform< TScalarType, NDimensions >
::GetFixedParameters(void) const
{
  unsigned int paramSize = 0, par = 0;

  for ( unsigned int t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    paramSize += m_LocalAffineTransformVector[t]->GetFixedParameters().GetSize();
    }

  this->m_FixedParameters.SetSize (paramSize);

  for ( unsigned int t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    const LocalAffineParametersType &param = m_LocalAffineTransformVector[t]->GetFixedParameters();
    for ( unsigned int par1 = 0; par1 < param.GetSize(); par1++ )
      {
      this->m_FixedParameters[par] = param[par1];
      par++;
      }
    }

  return this->m_FixedParameters;
}

// Get parameters
template< class TScalarType, unsigned int NDimensions >
const typename PolyAffineTransform< TScalarType,
                                          NDimensions >::ParametersType &
PolyAffineTransform< TScalarType, NDimensions >
::GetParameters(void) const
{
  unsigned int paramSize = 0, par = 0;

  for ( unsigned int t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    paramSize += m_LocalAffineTransformVector[t]->GetNumberOfParameters();
    }

  this->m_Parameters.SetSize (paramSize);

  for ( unsigned int t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    const LocalAffineParametersType &param = m_LocalAffineTransformVector[t]->GetParameters();
    for ( unsigned int par1 = 0; par1 < param.GetSize(); par1++ )
      {
      this->m_Parameters[par] = param[par1];
      par++;
      }
    }

  return this->m_Parameters;
}

// Set parameters
template< class TScalarType, unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::SetParameters(const ParametersType & parameters)
{
  unsigned int par = 0;

  //Save parameters. Needed for proper operation of TransformUpdateParameters.
  if( &parameters != &(this->m_Parameters) )
    {
    this->m_Parameters = parameters;
    }

  unsigned int paramSize;

  for ( unsigned int t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    //assuming each transform object has m_Parameters set initially
    LocalAffineParametersType param = m_LocalAffineTransformVector[t]->GetParameters();
    paramSize = param.GetSize();

    for ( unsigned int par1 = 0; par1 < paramSize; par1++ )
      {
      param[par1] = this->m_Parameters[par];
      par++;
      }

    m_LocalAffineTransformVector[t]->SetParameters(param);
    }

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();
}

// Compute the Jacobian in one position
template< class TScalarType,
          unsigned int NDimensions >
const typename PolyAffineTransform< TScalarType, NDimensions >::JacobianType &
PolyAffineTransform< TScalarType, NDimensions >
::GetJacobian(const InputPointType & p) const
{
  this->GetJacobianWithRespectToParameters( p, this->m_SharedLocalJacobian );
  return this->m_SharedLocalJacobian;

}

// Compute the Jacobian in one position, without setting values to m_Jacobian
template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::GetJacobianWithRespectToParameters(const InputPointType & p, JacobianType &j) const
{
  //This will not reallocate memory if the dimensions are equal
  // to the matrix's current dimensions.

  j.SetSize( NDimensions, this->GetNumberOfLocalParameters() );
  j.Fill(0.0);

  unsigned int t, d, par1, par = 0;
  double squaredDistance, *weights, sumWeight;

  sumWeight = 0;
  weights = new double[ m_LocalAffineTransformVector.size() ];

  for ( t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    squaredDistance = p.SquaredEuclideanDistanceTo(m_LocalAffineTransformVector[t]->GetCenter());
    weights[t] = vcl_exp( - squaredDistance );

    sumWeight = sumWeight + weights[t];
    }
  for ( t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    weights[t] /= sumWeight; //normalizing weights
    }

  for ( t = 0; t < m_LocalAffineTransformVector.size(); t++ )
    {
    const typename LocalAffineTransformType::JacobianType &atomJacobian
      = m_LocalAffineTransformVector[t]->GetJacobian(p);
    for ( par1 = 0; par1 < m_LocalAffineTransformVector[t]->GetNumberOfLocalParameters(); par1++ )
      {
      for ( d = 0; d < NDimensions; d++ )
        {
        j(d, par) = weights[t] * atomJacobian[d][par1];
        }
      par++;
      }
    }

  delete weights;

  return;

}

template< class TScalarType,
          unsigned int NDimensions >
TScalarType
PolyAffineTransform< TScalarType, NDimensions >
::GetDiagonalSpacing(MaskImagePointer mask)
{
  IndexType index0, index1;
  index0.Fill(0);
  index1.Fill(1);

  PointType point0, point1;
  mask->TransformIndexToPhysicalPoint(index0, point0);
  mask->TransformIndexToPhysicalPoint(index1, point1);

  return point1.EuclideanDistanceTo(point0);
}

/**
 * Compute the union of the image domains of all transforms.
 * It uses the first image domain's direction and origin, but
 * with a spacing from the first image scaled to the minimum
 * diagonal spacing among all transforms.
 *
 * For each transform, its moving image domain is used since
 * it includes and might be bigger than its fixed image domain.
 */
template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::InitializeBoundaryMask()
{
  ContinuousIndexType minIndex, maxIndex;
  PointType point;
  RegionType region;
  SpacingType firstSpacing;
  TScalarType diagSpacing, minDiagSpacing, firstDiagSpacing;
  MaskImagePointer mask, firstMask;

  //use the direction matrix from the first image domain
  firstMask = m_LocalAffineTransformVector[0]->GetMovingMaskImage();

  //compute the minimum diagonal spacing from all transforms
  firstSpacing = firstMask->GetSpacing();
  firstDiagSpacing = this->GetDiagonalSpacing(firstMask);
  minDiagSpacing = firstDiagSpacing;

  //compute minIndex and maxIndex of all mask image corners.
  //minIndex and maxIndex are in the first image domain.
  region = firstMask->GetLargestPossibleRegion();
  minIndex = region.GetIndex();
  maxIndex = region.GetUpperIndex();

  for (int t=0; t<m_LocalAffineTransformVector.size(); t++)
    {
    LocalAffineTransformPointer trans = m_LocalAffineTransformVector[t];
    mask = trans->GetMovingMaskImage();
  
    //the minimum spacing is computed here.
    diagSpacing = this->GetDiagonalSpacing(mask);
    if (minDiagSpacing > diagSpacing)
      {
      minDiagSpacing = diagSpacing;
      }

    region = mask->GetLargestPossibleRegion();
    typedef typename itk::VectorContainer<int, IndexType> IndexContainerType;
    typename IndexContainerType::Pointer corners = 
      PicslImageHelper::GetCorners<RegionType>(region);
    typename IndexContainerType::ConstIterator corner;
    for ( corner = corners->Begin(); corner != corners->End(); corner++)
      {
      IndexType cornerIndex = corner.Value();
      PointType cornerPoint;
      mask->TransformIndexToPhysicalPoint(cornerIndex, cornerPoint);

      //get corner index in the first mask image
      ContinuousIndexType mappedCorner;
      firstMask->TransformPhysicalPointToContinuousIndex(cornerPoint, mappedCorner);

      PicslImageHelper::CopyWithMin<ContinuousIndexType>(minIndex, mappedCorner);
      PicslImageHelper::CopyWithMax<ContinuousIndexType>(maxIndex, mappedCorner);
      }
    } //for each transform

  TScalarType scalingFactor = minDiagSpacing / firstDiagSpacing;
  SpacingType unionSpacing = firstSpacing * scalingFactor;

  ContinuousIndexType minScaledIndex, maxScaledIndex;
  for (int d=0; d<NDimensions; d++)
    {
    minScaledIndex[d] = minIndex[d] / scalingFactor;
    maxScaledIndex[d] = maxIndex[d] / scalingFactor;
    //pad with some boundaries
    minScaledIndex[d] -= this->m_PadBoundary;
    maxScaledIndex[d] += this->m_PadBoundary;
    }

  IndexType paddedIndex1, paddedIndex2;
  paddedIndex1.CopyWithRound(minScaledIndex);
  paddedIndex2.CopyWithRound(maxScaledIndex);

  RegionType unionRegion;
  unionRegion.SetIndex( paddedIndex1 );
  unionRegion.SetUpperIndex( paddedIndex2 );

  this->m_BoundaryMask = MaskImageType::New();
  this->m_BoundaryMask->SetOrigin(firstMask->GetOrigin());
  this->m_BoundaryMask->SetDirection(firstMask->GetDirection());
  this->m_BoundaryMask->SetSpacing( unionSpacing );
  this->m_BoundaryMask->SetRegions( unionRegion );

  this->m_BoundaryMask->Allocate();

  this->m_BoundaryMask->FillBuffer(1);

  //set its boundary to zeros
  typedef typename itk::NeighborhoodAlgorithm
    ::ImageBoundaryFacesCalculator<MaskImageType> FaceCalculatorType;
  typedef typename FaceCalculatorType::FaceListType FaceListType;
  FaceCalculatorType faceCalculator;
  FaceListType faceList;

  SizeType radius;
  radius.Fill(1);

  faceList = faceCalculator(this->m_BoundaryMask, unionRegion, radius);
  typename FaceListType::iterator fit = faceList.begin();
  ++fit; //skip the first non-boundary face
  for( ; fit != faceList.end(); ++fit)
    {
    itk::ImageRegionIterator<MaskImageType> it(this->m_BoundaryMask, *fit);
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
      it.Set(0); //set the boundary pixels to zero
      }
    }
}

template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::RewindTrajectory(unsigned int transformId, int stopTime)
{
  this->m_LocalAffineTransformVector[transformId]->SetStopTime(stopTime);

  TrajectoryImagePointer traj = this->m_TrajectoryImageVector[transformId];
  int stopValue = this->TimeStampToImageValue(stopTime);

  typedef ImageRegionIteratorWithIndex< TrajectoryImageType > IteratorType;
  IteratorType it( traj, traj->GetLargestPossibleRegion() );

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    if (it.Get() > stopValue)
      {
      it.Set(0); //set to background
      }
    }  
}

/**
 * Compute the trajectory of the spatial object under each transform by
 * propagation at each time stamp.
 */
template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::ComputeNextStepTrajectory(unsigned int transformId)
{
  this->m_TimerTrajectory.Start();

  bool overlap = false;
  unsigned int overlapPointId = 0;
  LocalAffineTransformPointer trans = this->m_LocalAffineTransformVector[transformId];

  int timeStamp = trans->GetStopTime();
  if (timeStamp >= trans->GetTimeStampMax())
    {
    this->m_TimerTrajectory.Stop();
    return;
    }

  timeStamp++;
  trans->SetStopTime(timeStamp);

  typedef typename LocalAffineTransformType::AffineTransformPointer AffineTransformPointer;
  AffineTransformPointer partialTrans = trans->GetPartialTransform(0, timeStamp);

  TrajectoryImagePointer traj = this->m_TrajectoryImageVector[transformId];
  PointSetPointer samples = trans->GetSamplePointSet();
  unsigned int pointId = 0;

  for (int f=0; f<samples->GetNumberOfPoints(); f++)
    {
		PointType z, x = samples->GetPoint(f);
    IndexType index;

		//compute z = exp(t*log(T) * x
    z = partialTrans->TransformPoint( x );
    bool insideImage = traj->TransformPhysicalPointToIndex(z, index);

    if (insideImage)
      {
      if (traj->GetPixel(index) == 0) //not in trajectory yet
        {
        traj->SetPixel(index, this->TimeStampToImageValue(timeStamp));

        //overlap detection
        overlap = PointExistsInOtherTrajectories(transformId, index);
        if (overlap)
          {
            overlapPointId = f; //not used now, may be used later
            break;
          }
        }
      } //if insideImage
    } //end samples

  trans->SetOverlapped(overlap);
  trans->SetOverlapPointId(overlapPointId);

  this->m_TimerTrajectory.Stop();
}

template< class TScalarType,
          unsigned int NDimensions >
bool
PolyAffineTransform< TScalarType, NDimensions >
::PointExistsInOtherTrajectories(unsigned int transformId, IndexType index)
{
  TrajectoryImagePointer trajOther;
  bool exist = false;

  for (unsigned int t=0; t<this->GetNumberOfLocalAffineTransforms(); t++)
    {
    if (t == transformId)
      {
      continue;
      }
    trajOther = this->m_TrajectoryImageVector[t];
    if (trajOther->GetPixel(index) != 0) //already in another trajectory
      {
      exist = true;
      break;
      }
    } //end for other trajectories

  return exist;
}

template< class TScalarType,
          unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::DistanceMapImagePointer
PolyAffineTransform< TScalarType, NDimensions >
::ComputeTrajectoryDistanceMapImage(TrajectoryImagePointer traj)
{
  this->m_TimerDistanceMap.Start();

  typedef itk::SignedMaurerDistanceMapImageFilter
    <TrajectoryImageType, DistanceMapImageType>                DistanceMapImageFilterType;
  typedef typename DistanceMapImageFilterType::Pointer         DistanceMapImageFilterPointer;

  DistanceMapImageFilterPointer filter = DistanceMapImageFilterType::New();
    
  filter->SetInput( traj );
  filter->SetSquaredDistance( false );
  filter->SetUseImageSpacing( true );
  filter->SetInsideIsPositive( false ); //outside is positive
  filter->Update();

  DistanceMapImagePointer output = filter->GetOutput();

  this->m_TimerDistanceMap.Stop();

  return output;
}

template< class TScalarType,
          unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::DistanceMapImagePointer
PolyAffineTransform< TScalarType, NDimensions >
::ComputeBoundaryDistanceMapImage()
{
  typedef itk::SignedMaurerDistanceMapImageFilter
    <MaskImageType, DistanceMapImageType>                           DistanceMapImageFilterType;
  typedef typename DistanceMapImageFilterType::Pointer         DistanceMapImageFilterPointer;

  DistanceMapImageFilterPointer filter = DistanceMapImageFilterType::New();
    
  filter->SetInput( this->m_BoundaryMask );
  filter->SetSquaredDistance( false );
  filter->SetUseImageSpacing( true );
  filter->SetInsideIsPositive( true ); //inside is positive
  filter->Update();

  DistanceMapImagePointer output = filter->GetOutput();
  return output; 
}

/** This method initializes the buffers. 
 *  
 *  (1) For each local affine transform, it dilates the fixed mask, 
 *  compute the moving mask produced by the transform, and compute
 *  the sample point set at time 0.
 *
 *  (2) It compute m_BoundaryMask as the union of the image
 *  masks of all transforms. Then it computes the distance map
 *  of the boundary mask.
 *
 *  (3) It allocates the array for trajectories of local transforms
 *  and the array for distance maps of the trajectories. Also it 
 *  allocates the memory buffer for each trajectory and the combined
 *  trajectory. Note: the memory buffers for distance maps are 
 *  allocated in the distance map filter.
 *
 *  (4) It allocates the memory buffers for m_DisplacementField
 *  and m_VelocityField.
 */
template< class TScalarType,
          unsigned int NDimensions >
bool
PolyAffineTransform< TScalarType, NDimensions >
::InitializeBuffers()
{
  m_TimerInitializeBuffers.Start();

  if (this->GetNumberOfLocalAffineTransforms() == 0)
    {
    m_TimerInitializeBuffers.Stop();
    itkExceptionMacro(
      << "The number of local affine transforms is zero." );
    return false;
    }
  
  for (unsigned int t=0; t<this->GetNumberOfLocalAffineTransforms(); t++)
    {
    LocalAffineTransformPointer trans = this->m_LocalAffineTransformVector[t];

    trans->DilateFixedMaskImage(this->m_PadTrajectory);
    trans->ComputeMovingMaskImage();  
    trans->ComputeSamplePointSet(0);

    PicslImageHelper::WriteImage<MaskImageType>(
      trans->GetFixedMaskImage(), "tmpFixedMask.nii", (int)(t+1));
    PicslImageHelper::WriteImage<MaskImageType>(
      trans->GetMovingMaskImage(), "tmpMovingMask.nii", (int)(t+1));
    }

  //Initialize BoundaryMask and its DistanceMapImage
  this->InitializeBoundaryMask();
  PicslImageHelper::WriteImage<MaskImageType>
    (this->m_BoundaryMask, "tmpBoundMask.nii");

  this->m_BoundaryDistanceMapImage = this->ComputeBoundaryDistanceMapImage();
  PicslImageHelper::WriteImage<DistanceMapImageType>
    (this->m_BoundaryDistanceMapImage, "tmpBoundDistance.nii");

  this->UpdateTimeStampIfAutomatic();
  for (unsigned int t=0; t<this->GetNumberOfLocalAffineTransforms(); t++)
    {
    LocalAffineTransformPointer trans = this->m_LocalAffineTransformVector[t];
    trans->SetTimeStampMax(this->m_TimeStampMax);
    }

  //Initialize trajectory vector, its elments, and its frontier point sets
  unsigned int transformNumber = this->GetNumberOfLocalAffineTransforms();
  this->m_TrajectoryImageVector.resize(transformNumber);
  this->m_TrajectoryDistanceMapImageVector.resize(transformNumber);

  TrajectoryImagePointer traj = TrajectoryImageType::New();
  traj->CopyInformation(this->m_BoundaryMask);
  traj->SetRegions(this->m_BoundaryMask->GetLargestPossibleRegion());
  traj->Allocate();
  this->m_CombinedTrajectoryImage = traj;

  for (unsigned int t=0; t<this->GetNumberOfLocalAffineTransforms(); t++)
    {
    traj = TrajectoryImageType::New();
    traj->CopyInformation(this->m_BoundaryMask);
    traj->SetRegions(this->m_BoundaryMask->GetLargestPossibleRegion());
    traj->Allocate();
    this->m_TrajectoryImageVector[t] = traj;
    }

  //Initialize the overall displacement field
  this->m_DisplacementField = DisplacementFieldType::New();
  this->m_DisplacementField->CopyInformation(this->m_BoundaryMask);
  this->m_DisplacementField->SetRegions(this->m_BoundaryMask->GetLargestPossibleRegion());
  this->m_DisplacementField->Allocate();

  typedef typename DisplacementFieldType::PixelType VectorType;
  VectorType disp(0.0);
  this->m_DisplacementField->FillBuffer(disp); //identity transform

  //Initialize the sum of velocity field
  this->m_VelocityField = DisplacementFieldType::New();
  this->m_VelocityField->CopyInformation(this->m_BoundaryMask);
  this->m_VelocityField->SetRegions(this->m_BoundaryMask->GetLargestPossibleRegion());
  this->m_VelocityField->Allocate();

  m_TimerInitializeBuffers.Stop();
  return true;
}

template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::ComputeWeightedSumOfVelocityFields()
{
  IndexType index;
  PointType point;
  DisplacementVectorType velocitySum;

  unsigned int numTrans = this->GetNumberOfLocalAffineTransforms();
  int ownerTraj;
  double *distances = new double[numTrans];
  double distCombinedTraj, distBoundary;

  for (unsigned int t=0; t<numTrans; t++)
    {
    this->m_LocalAffineTransformVector[t]->UpdatePartialVelocityMatrix();
    }

  typedef ImageRegionIteratorWithIndex< DisplacementFieldType > IteratorType;
  IteratorType it( this->m_VelocityField, this->m_VelocityField->GetLargestPossibleRegion() );

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    this->m_VelocityField->TransformIndexToPhysicalPoint( index, point );

    ownerTraj = -1;
    for (unsigned int t=0; t<numTrans; t++)
      {
      distances[t] = this->m_TrajectoryDistanceMapImageVector[t]->GetPixel(index);	
      if (distances[t] <= 0)
        {
        distances[t] = 0;
        ownerTraj = t;
        }
      }
    if (index[0]==151-81 && index[1]==191-91)
      int debug=1;
    if (ownerTraj >= 0) //the point is inside a trajectory ownerTraj
      {
      LocalAffineTransformPointer localAffine = this->m_LocalAffineTransformVector[ownerTraj];
      velocitySum = localAffine->GetPartialVelocityAtPoint(point);
      }
    else //the point is outside all trajectories
      {
      distCombinedTraj = this->m_CombinedTrajectoryDistanceMapImage->GetPixel(index); 
      distBoundary = this->m_BoundaryDistanceMapImage->GetPixel(index);
      this->ComputeWeightedSumOfVelocitiesAtPoint(velocitySum, point,
        distances, distCombinedTraj, distBoundary);
      }

    it.Set(velocitySum);
    } //for it

  delete distances;
}

template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::ComputeWeightedSumOfVelocitiesAtPoint(DisplacementVectorType &velocitySum, const PointType &point,
  const double distances[], double distanceToCombinedTraj, double distanceToImageBoundary)
{
  unsigned int numTrans = this->GetNumberOfLocalAffineTransforms();

  //eliminator the common factor to avoid a very small denominator
  double distMin = NumericTraits<double>::max();
  for (unsigned int t=0; t<numTrans; t++)
    {
    if (distMin > distances[t])
      {
      distMin = distances[t];
      }
    }

  double weightSum = 0.0, weight;
  velocitySum.Fill(0.0);
  
  for (unsigned int t=0; t<numTrans; t++)
    {
    LocalAffineTransformPointer localAffine = this->m_LocalAffineTransformVector[t];
    DisplacementVectorType velocity = localAffine->GetPartialVelocityAtPoint(point);

    //weight = vcl_exp(- (distances[t]*distances[t]-distMin*distMin) / this->m_DecayRateOfTrajectory);
    weight = vcl_exp(- (distances[t]-distMin) / this->m_DecayRateOfTrajectory);
    weightSum += weight;
    velocitySum += velocity * weight;
    }

  if (distanceToCombinedTraj < 0)
    {
    distanceToCombinedTraj = 0;
    }
  if (distanceToImageBoundary < 0)
    {
    distanceToImageBoundary = 0;
    }
  double weightCombinedTraj = vcl_exp(
    //- distanceToCombinedTraj*distanceToCombinedTraj / this->m_DecayRateOfTrajectory);
    - distanceToCombinedTraj / this->m_DecayRateOfTrajectory);
  double weightBoundary = 1 - vcl_exp(
    //- distanceToImageBoundary*distanceToImageBoundary / this->m_DecayRateOfBoundary);
    - distanceToImageBoundary / this->m_DecayRateOfBoundary);

  velocitySum *= (1.0 / weightSum * weightCombinedTraj * weightBoundary);
}

/** Compute the minimum stop time stamps of all local transforms. 
 *  It also prints the stop time stamps if debugInfo is not null.
 */
template< class TScalarType,
          unsigned int NDimensions >
int
PolyAffineTransform< TScalarType, NDimensions >
::GetMinStopTime(char *debugInfo)
{
  int minStopTime = NumericTraits<int>::max();
  int stopTime;
  if (debugInfo)
    {
    std::cout << "GetMinStopTime: " << debugInfo;
    }
  for (unsigned int t = 0; t < this->GetNumberOfLocalAffineTransforms(); t++)
    {
    stopTime = this->m_LocalAffineTransformVector[t]->GetStopTime();
    if (debugInfo)
      {
      std::cout << " " << stopTime;
      }
    if (minStopTime > stopTime)
      {
      minStopTime = stopTime;
      }
    }
  if (debugInfo)
    {
    std::cout << std::endl;
    }
  return minStopTime;
}

/** Compute the velocity field before the trajectories overlap.
 *
 *  First it computes the trajectories of each local transform
 *  from its current m_StartTime until the overlap happens.
 *  There is a switch m_StopAllTrajectoriesAtOverlap which
 *  decides the behaviour when an overlap happens. If the swith
 *  is false, other trajectories without overlap will continue
 *  to move forward. Otherwise, all trajectories will stop.
 *
 *  Second it will rewind backwards the trajectories with overlap.
 *  Therefore, there will be some extra space between trajectories.
 *
 *  Second it compute the distance maps of these trajectories.
 *
 *  Third it compute the weighted sum of the velocity fields
 *  according to the distances.
 */
template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::ComputeVelocityFieldBeforeOverlap()
{
  m_TimerComputeVelocityFieldBeforeOverlap.Start();
  for (unsigned int t = 0; t < this->GetNumberOfLocalAffineTransforms(); t++)
    {
    this->m_LocalAffineTransformVector[t]->SetOverlapped(false);
    }

  bool stopNow = false;
  while (!stopNow && this->GetMinStopTime("BeforeOverlap") < this->m_TimeStampMax)
    {
    bool moved = false;
    for (unsigned int t = 0; !stopNow && t < this->GetNumberOfLocalAffineTransforms(); t++)
      {
      LocalAffineTransformPointer trans = this->m_LocalAffineTransformVector[t];
      if (!trans->GetOverlapped() && trans->GetStopTime() < this->m_TimeStampMax)
        {
        this->ComputeNextStepTrajectory(t);
        moved = true;
        }
      if (trans->GetOverlapped() && this->m_StopAllTrajectoriesAtOverlap)
        {
        stopNow = true;
        }
      }
    if (!moved)
      {
      break;
      }
    }

  // save traj for debug
  for (unsigned int t = 0; t < this->GetNumberOfLocalAffineTransforms(); t++)
    {
    LocalAffineTransformPointer trans = this->m_LocalAffineTransformVector[t];

    PicslImageHelper::WriteImage<TrajectoryImageType>
      (this->m_TrajectoryImageVector[t], "tmpTraj.nii", 
        (t+1)*10000 + trans->GetStopTime());

    if (trans->GetOverlapped()) //rewind trajectory
      {
      int rewindTime = (int)(0.3*trans->GetStartTime() + 0.7*trans->GetStopTime()); //any better solution?
      this->RewindTrajectory(t, rewindTime); //rewind trajectory

      PicslImageHelper::WriteImage<TrajectoryImageType>
        (this->m_TrajectoryImageVector[t], "tmpTrajRewind.nii", 
        (t+1)*10000+rewindTime);
      }
    }

  for (unsigned int t = 0; t < this->GetNumberOfLocalAffineTransforms(); t++)
    {
    this->m_TrajectoryDistanceMapImageVector[t] = this->ComputeTrajectoryDistanceMapImage(
      this->m_TrajectoryImageVector[t]);

    PicslImageHelper::WriteImage<DistanceMapImageType>
      (this->m_TrajectoryDistanceMapImageVector[t], "tmpTrajDistance.nii", 
      (t+1)*10000 + this->m_LocalAffineTransformVector[t]->GetStopTime());
    }

  this->CombineTrajectories();
  this->m_CombinedTrajectoryDistanceMapImage = this->ComputeTrajectoryDistanceMapImage(
      this->m_CombinedTrajectoryImage);

  PicslImageHelper::WriteImage<TrajectoryImageType>
        (this->m_CombinedTrajectoryImage, "tmpTrajComb.nii", this->GetMinStopTime());  
  PicslImageHelper::WriteImage<DistanceMapImageType>
        (this->m_CombinedTrajectoryDistanceMapImage, "tmpTrajCombDistance.nii", this->GetMinStopTime());

  this->ComputeWeightedSumOfVelocityFields();

  m_TimerComputeVelocityFieldBeforeOverlap.Stop();
}

template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::CombineTrajectories()
{
  this->m_CombinedTrajectoryImage->FillBuffer(0);

  typedef ImageRegionIteratorWithIndex< TrajectoryImageType > IteratorType;

  for (unsigned int t=0; t<this->GetNumberOfLocalAffineTransforms(); t++)
    {
    TrajectoryImagePointer traj = this->m_TrajectoryImageVector[t];
    IteratorType it( traj, traj->GetLargestPossibleRegion() );

    for( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      int ts = it.Get();
      if (ts != 0)
        {
        this->m_CombinedTrajectoryImage->SetPixel(it.GetIndex(), ts);
        }
      }
    }
}

/** Get the displacement field of this PolyAffineTransform.
 */
template< class TScalarType,
          unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::DisplacementFieldType*
PolyAffineTransform< TScalarType, NDimensions >
::GetDisplacementField()
{
  if (this->m_DisplacementField.IsNull())
    {
    this->ComputeDisplacementField();
    }
  
  this->PrintTimers();

  return this->m_DisplacementField;
}

/** Compute the displacement field of this PolyAffineTransform.
 *  (1) It initializes the trajectories of local transforms.
 *  (2) It will compute the velocity sum for each segment of time
 *  during which there is no trajectory overlap. 
 *  (3) It computes the expontial mapping to get a displacement field
 *  during that time segment.
 *  (4) It computes the composition of these displacement fields at
 *  mutiple tiem segments.
 */
template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::ComputeDisplacementField()
{
  this->m_TimerComputeDisplacementField.Start();

  this->InitializeBuffers();

  for (unsigned int t=0; t<this->GetNumberOfLocalAffineTransforms(); t++)
    {
    this->m_LocalAffineTransformVector[t]->SetStopTime(0);
    this->m_LocalAffineTransformVector[t]->SetTimeStampMax(this->m_TimeStampMax);
    }

  while (this->GetMinStopTime("MainLoop") < this->m_TimeStampMax)
    {
    m_TimerInitializeIteration.Start();
    for (unsigned int t=0; t<this->GetNumberOfLocalAffineTransforms(); t++)
      {
      LocalAffineTransformPointer trans = this->m_LocalAffineTransformVector[t];

      // The next trajectory always starts at the m_StopTime of the previous trajectory.
      int DonePreviously = trans->GetStopTime();
      trans->SetStartTime(DonePreviously);

      // For each local transform, its trajectory is during [ m_StartTime, m_StopTime ].
      // at the beginning, each trajectory is empty and m_StopTime = m_StartTime - 1.
      trans->SetStopTime(DonePreviously - 1);

      this->m_TrajectoryImageVector[t]->FillBuffer(0);
      }
    m_TimerInitializeIteration.Stop();

    //Compute the velocity field before overlap. It computes the trajectories,
    //related weights, and weighted sum of velocities.
    this->ComputeVelocityFieldBeforeOverlap();

    DisplacementFieldPointer exponentialField =
      this->ComputeExponentialDisplacementField(this->m_VelocityField);

    PicslImageHelper::WriteImage<DisplacementFieldType>
      (this->m_VelocityField, "tmpVelocitySum.nii", this->GetMinStopTime());
    PicslImageHelper::WriteDisplacementField<DisplacementFieldType>
      (this->m_VelocityField, "tmpVelocitySum.nii", this->GetMinStopTime());
    PicslImageHelper::WriteDisplacementField<DisplacementFieldType>
      (exponentialField, "tmpExponentField.nii", this->GetMinStopTime());
 
    typedef itk::ComposeDisplacementFieldsImageFilter<DisplacementFieldType> ComposerType;
    typename ComposerType::Pointer composer = ComposerType::New();
    composer->SetDisplacementField( exponentialField );
    composer->SetWarpingField( this->m_DisplacementField );
    this->m_TimerDisplacementFieldComposing.Start();
    composer->Update();
    this->m_TimerDisplacementFieldComposing.Stop();

    this->m_DisplacementField = composer->GetOutput();
    this->m_DisplacementField->DisconnectPipeline();
    }

  this->m_TimerComputeDisplacementField.Stop();
}

template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::UpdateTimeStampIfAutomatic()
{
  if (this->m_BoundaryMask.IsNull())
    {
    itkExceptionMacro("Please initialize m_BoundaryMask before automatically setting m_TimeStampMax.");
    }

  if (this->m_TimeStampAutomaticFlag)
    {
    RegionType region = this->m_BoundaryMask->GetLargestPossibleRegion();
    SizeType size = region.GetSize();

    double diagonalLength = 0.0;
    for (unsigned int d=0; d<NDimensions; d++)
      {
      diagonalLength += size[d] * size[d];
      }
    diagonalLength = vcl_sqrt(diagonalLength);

    int diagonalVoxels = itk::Math::RoundHalfIntegerUp<int, double>(diagonalLength);

    int nonzeroBit, leftmostBit = -1;
    for (int bit=0; bit<sizeof(diagonalVoxels)*8; bit++)
      {
      nonzeroBit = 1 << bit;
      if ( (diagonalVoxels & nonzeroBit) != 0)
        {
        leftmostBit = bit;
        }
      }

    this->m_TimeStampLog = leftmostBit + 1;
    std::cout << "Automatic m_TimeStampLog = " << this->m_TimeStampLog << std::endl;

    //m_TimeStampMax is the number of time stamps
    this->m_TimeStampMax = 1 << this->m_TimeStampLog;
    }
}


template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::PrintTimers()
{
  double totalTimerComputeFixedMaskImage = 0;
  double totalTimerComputeMovingMaskImage = 0;
  double totalTimerMatrixExponential = 0;
  double totalTimerMatrixLogarithm = 0;

  for (unsigned int t=0; t<this->GetNumberOfLocalAffineTransforms(); t++)
    {
    LocalAffineTransformPointer trans = this->m_LocalAffineTransformVector[t];
    totalTimerComputeFixedMaskImage += trans->GetTimerComputeFixedMaskImage().GetTotal();
    totalTimerComputeMovingMaskImage += trans->GetTimerComputeMovingMaskImage().GetTotal();
    totalTimerMatrixExponential += trans->GetTimerMatrixExponential().GetTotal();
    totalTimerMatrixLogarithm += trans->GetTimerMatrixLogarithm().GetTotal();
    }
  std::cout << "totalTimerComputeFixedMaskImage = " 
    << totalTimerComputeFixedMaskImage << " seconds" << std::endl;
  std::cout << "totalTimerComputeMovingMaskImage = " 
    << totalTimerComputeMovingMaskImage << " seconds" << std::endl;

  std::cout << "totalTimerMatrixExponential = " 
    << totalTimerMatrixExponential << " seconds" << std::endl;
  std::cout << "totalTimerMatrixLogarithm = " 
    << totalTimerMatrixLogarithm << " seconds" << std::endl;

  std::cout << std::endl;
  std::cout << "m_TimerTrajectory.GetTotal() = " 
    << m_TimerTrajectory.GetTotal() << " seconds" << std::endl;
  std::cout << "m_TimerDistanceMap.GetTotal() = " 
    << m_TimerDistanceMap.GetTotal() << " seconds" << std::endl;

  std::cout << "m_TimerExponentialMapping.GetTotal() = " 
    << m_TimerExponentialMapping.GetTotal() << " seconds" << std::endl;
  std::cout << "m_TimerDisplacementFieldComposing.GetTotal() = " 
    << m_TimerDisplacementFieldComposing.GetTotal() << " seconds" << std::endl;

  std::cout << std::endl;
  std::cout << "m_TimerComputeVelocityFieldBeforeOverlap.GetTotal() = " 
    << m_TimerComputeVelocityFieldBeforeOverlap.GetTotal() << " seconds" << std::endl;
  std::cout << "m_TimerInitializeBuffers.GetTotal() = " 
    << m_TimerInitializeBuffers.GetTotal() << " seconds" << std::endl;
  std::cout << "m_TimerInitializeIteration.GetTotal() = " 
    << m_TimerInitializeIteration.GetTotal() << " seconds" << std::endl;

  std::cout << std::endl;
  std::cout << "m_TimerComputeDisplacementField.GetTotal() = " 
    << m_TimerComputeDisplacementField.GetTotal() << " seconds" << std::endl;
  
}
template< class TScalarType,
          unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::DisplacementFieldPointer
PolyAffineTransform< TScalarType, NDimensions >
::ComputeExponentialDisplacementField(const DisplacementFieldPointer &velocityField)
{
  this->m_TimerExponentialMapping.Start();

  ExponentialImageFilterPointer filter = ExponentialImageFilterType::New();
  filter->SetInput(velocityField);

  filter->SetAutomaticNumberOfIterations(true);
  filter->Update();
  DisplacementFieldPointer exponentialField = filter->GetOutput();
  
  this->m_TimerExponentialMapping.Stop();

  return exponentialField;
}

} // namespace

#endif
