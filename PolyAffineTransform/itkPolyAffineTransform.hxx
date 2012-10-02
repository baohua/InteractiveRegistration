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
  this->m_DecayConstant = 1.0/8;
  this->SetTimeStampLog(8); //by default 2^8 time stamps
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
    minScaledIndex[d]--;
    maxScaledIndex[d]++;
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

/**
 * Initialize the frontier during trajectory computation. The frontier
 * contains the points that are in the trajectory, but not have been
 * processed for its move at the next time stamp.
 */
template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::InitializeFrontier(unsigned int transformId)
{
  PointSetPointer pointSet = m_LocalAffineTransformVector[transformId]->GetDenseFixedPointSet();
  SizeValueType num = pointSet->GetNumberOfPoints();

  PointSetPointer frontier = PointSetType::New(); 

  typename PointSetType::PointsContainer *container = pointSet->GetPoints();
  typename PointSetType::PointsContainer::ConstIterator it = container->Begin();

  unsigned int pointId = 0;
  while (it != container->End())
    {
    frontier->SetPoint(pointId++, it->Value());
    it++;
    }

  m_FrontierVector[transformId] = frontier;
}

template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::InitializeTrajectory(unsigned int transformId)
{
  TrajectoryImagePointer traj = TrajectoryImageType::New();
  traj->CopyInformation(this->m_BoundaryMask);
  traj->SetRegions(this->m_BoundaryMask->GetLargestPossibleRegion());
  traj->Allocate();
  traj->FillBuffer(0);

  this->m_TrajectoryImageVector[transformId] = traj;
}

template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::RewindTrajectory(unsigned int transformId, int stopTime)
{
  TrajectoryImagePointer traj = this->m_TrajectoryImageVector[transformId];
  int stopTimeStamp = this->TranslateTimeStamp(stopTime);

  typedef ImageRegionIteratorWithIndex< TrajectoryImageType > IteratorType;
  IteratorType it( traj, traj->GetLargestPossibleRegion() );

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    if (it.Get() > stopTimeStamp)
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
::ComputeNextTimeTrajectory(unsigned int transformId, int timeStamp,
                            bool &overlap, unsigned int &overlapPointId)
{
  LocalAffineTransformPointer trans = this->m_LocalAffineTransformVector[transformId];

  typedef typename LocalAffineTransformType::VelocityAffineTransformPointer VelocityAffineTransformPointer;
  VelocityAffineTransformPointer velocity = trans->GetVelocityAffineTransform();

  TrajectoryImagePointer traj = this->m_TrajectoryImageVector[transformId];
  TrajectoryImagePointer trajOther;

  PointSetPointer frontierNow = this->m_FrontierVector[transformId];
  PointSetPointer frontierNext = PointSetType::New();
  unsigned int pointId = 0;

  overlap = false;
  for (int f=0; !overlap && f<frontierNow->GetNumberOfPoints(); f++)
    {
		PointType y = frontierNow->GetPoint(f);
    PointType z, step;

		//compute z=(1+V/N)*y
    step = velocity->TransformPoint( y );
    for (unsigned int d=0; d<NDimensions; d++)
      {
      z[d] = y[d] + step[d] / this->m_TimeStampExp;
      }
    frontierNext->SetPoint(pointId++, z);

    IndexType index;
    bool insideImage = traj->TransformPhysicalPointToIndex(y, index);
    if (insideImage)
      {
      if (traj->GetPixel(index) == 0)
        {
        traj->SetPixel(index, this->TranslateTimeStamp(timeStamp));

        //overlap detection
        for (unsigned int t=0; !overlap && t<this->m_LocalAffineTransformVector.size(); t++)
          {
          if (t == transformId)
            {
            continue;
            }
          trajOther = this->m_TrajectoryImageVector[t];
          if (trajOther->GetPixel(index) != 0) //already in another trajectory
            {
            overlap = true;
            overlapPointId = f;
            }
          } //end for other trajectories
        }
      }
    }
  
  this->m_FrontierVector[transformId] = frontierNext;

}

template< class TScalarType,
          unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::WeightImagePointer
PolyAffineTransform< TScalarType, NDimensions >
::ComputeTrajectoryWeightImage(TrajectoryImagePointer traj,
                               WeightImagePointer boundaryWeightImage)
{
  typedef itk::SignedMaurerDistanceMapImageFilter
    <TrajectoryImageType, WeightImageType>                     DistanceMapImageFilterType;
  typedef typename DistanceMapImageFilterType::Pointer         DistanceMapImageFilterPointer;

  DistanceMapImageFilterPointer filter = DistanceMapImageFilterType::New();
    
  filter->SetInput( traj );
  filter->SetSquaredDistance( false );
  filter->SetUseImageSpacing( true );
  filter->SetInsideIsPositive( false );
  filter->Update();

  WeightImagePointer output = filter->GetOutput();
  TScalarType weight, distance, boundaryWeight;
  IndexType index;

  typedef ImageRegionIteratorWithIndex< WeightImageType > IteratorType;
  IteratorType it( output, output->GetLargestPossibleRegion() );

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    distance = it.Get();
    index = it.GetIndex();
    boundaryWeight = boundaryWeightImage->GetPixel(index);

    if (distance > 0) 
      {
      weight = vcl_exp(-distance * this->m_DecayConstant); //weight decreases with distance
      weight *= boundaryWeight;
      }
    else
      {
      weight = 1;
      }
    it.Set(weight);
    }

  return output;
}

template< class TScalarType,
          unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::WeightImagePointer
PolyAffineTransform< TScalarType, NDimensions >
::ComputeBoundaryWeightImage()
{
  typedef itk::SignedMaurerDistanceMapImageFilter
    <MaskImageType, WeightImageType>                           DistanceMapImageFilterType;
  typedef typename DistanceMapImageFilterType::Pointer         DistanceMapImageFilterPointer;

  DistanceMapImageFilterPointer filter = DistanceMapImageFilterType::New();
    
  filter->SetInput( this->m_BoundaryMask );
  filter->SetSquaredDistance( false );
  filter->SetUseImageSpacing( true );
  filter->SetInsideIsPositive( false );
  filter->Update();

  WeightImagePointer output = filter->GetOutput();
  TScalarType weight, distance;

  typedef ImageRegionIteratorWithIndex< WeightImageType > IteratorType;
  IteratorType it( output, output->GetLargestPossibleRegion() );

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    distance = it.Get();

    if (distance < 0) //inside points have negative distance
      {
      //std::cout << "boundary distance = " << distance << std::endl;
      weight = 1 - vcl_exp(distance * this->m_DecayConstant); //weight increases with abs(distance)
      }
    else
      {
      weight = 0;
      }
    it.Set(weight);
    }

  return output;
}

template< class TScalarType,
          unsigned int NDimensions >
bool
PolyAffineTransform< TScalarType, NDimensions >
::InitializeBuffers()
{
  if (this->m_LocalAffineTransformVector.size() == 0)
    {
    itkWarningMacro(
      << "The number of local affine transforms is zero." );
    return false;
    }

  //Initialize BoundaryMask and its WeightImage
  this->InitializeBoundaryMask();
  PicslImageHelper::WriteImage<MaskImageType>
    (this->m_BoundaryMask, "tmpBoundMask.nii");

  this->m_BoundaryWeightImage = this->ComputeBoundaryWeightImage();
  PicslImageHelper::WriteImage<WeightImageType>
    (this->m_BoundaryWeightImage, "tmpBoundWeight.nii");

  //Initialize trajectory vector, its elments, and its frontier point sets
  unsigned int transformNumber = this->m_LocalAffineTransformVector.size();
  this->m_TrajectoryWeightImageVector.resize(transformNumber);
  this->m_TrajectoryImageVector.resize(transformNumber);
  this->m_FrontierVector.resize(transformNumber);

  for (unsigned int t=0; t<transformNumber; t++)
    {
    this->InitializeTrajectory(t);
    this->InitializeFrontier(t);
    }

  //Initialize the overall displacement field
  this->m_DisplacementField = DisplacementFieldType::New();
  this->m_DisplacementField->CopyInformation(this->m_BoundaryMask);
  this->m_DisplacementField->SetRegions(this->m_BoundaryMask->GetLargestPossibleRegion());
  this->m_DisplacementField->Allocate();

  //Initialize the sum of velocity field
  this->m_VelocityField = DisplacementFieldType::New();
  this->m_VelocityField->CopyInformation(this->m_BoundaryMask);
  this->m_VelocityField->SetRegions(this->m_BoundaryMask->GetLargestPossibleRegion());
  this->m_VelocityField->Allocate();

  return true;
}

template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::ComputeWeightedSumOfVelocityFields()
{
  typedef typename LocalAffineTransformType::VelocityAffineTransformType 
    VelocityAffineTransformType;
  typedef typename LocalAffineTransformType::VelocityAffineTransformPointer
    VelocityAffineTransformPointer;

  IndexType index;
  PointType point, velocity;
  TScalarType weight;
  typename DisplacementFieldType::PixelType velocitySum;

  typedef ImageRegionIteratorWithIndex< DisplacementFieldType > IteratorType;
  IteratorType it( this->m_VelocityField, this->m_VelocityField->GetLargestPossibleRegion() );

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    this->m_VelocityField->TransformIndexToPhysicalPoint( index, point );

    velocitySum.Fill(0.0);
    for (int t=0; t<this->m_LocalAffineTransformVector.size(); t++)
      {
      LocalAffineTransformPointer localAffine = this->m_LocalAffineTransformVector[t];
      VelocityAffineTransformPointer localVelocity = localAffine->GetVelocityAffineTransform();

      velocity = localVelocity->TransformPoint(point);
      weight = this->m_TrajectoryWeightImageVector[t]->GetPixel(index);

      for (unsigned int j=0; j<NDimensions; j++)
        {
        velocitySum[j] += weight * velocity[j];
        }
      }
    it.Set(velocitySum);
    }
}

template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::ComputeVelocityFieldBeforeOverlap(int &stopTime)
{
  bool overlap = false;
  int overlapTime;
  unsigned int overlapPointId, overlapTransformId;

  for (int ts = this->m_TimeStampMin; !overlap && ts <= this->m_TimeStampMax; ts++)
    {
    for (unsigned int t = 0; !overlap && t < this->m_LocalAffineTransformVector.size(); t++)
      {
      this->ComputeNextTimeTrajectory(t, ts, overlap, overlapPointId);
      if (overlap)
        {
        overlapTime = ts;
        overlapTransformId = t;
        }
      }
    }

  if (!overlap)
    {
    overlapTime = this->m_TimeStampMax;
    stopTime = this->m_TimeStampMax;
    }
  for (unsigned int t = 0; t < this->m_LocalAffineTransformVector.size(); t++)
    {
    PicslImageHelper::WriteImage<TrajectoryImageType>
      (this->m_TrajectoryImageVector[t], "tmpTraj.nii", t*10000+overlapTime);
    }

  if (overlap) //rewind trajectory
    {
    stopTime = overlapTime - 10; //any better solution?
    if (stopTime < this->m_TimeStampMin)
      {
      stopTime = this->m_TimeStampMin;
      }
    for (unsigned int t = 0; t < this->m_LocalAffineTransformVector.size(); t++)
      {
      this->RewindTrajectory(t, stopTime);
      PicslImageHelper::WriteImage<TrajectoryImageType>
        (this->m_TrajectoryImageVector[t], "tmpTrajRewind.nii", t*10000+stopTime);
      }
    }

  for (unsigned int t = 0; t < this->m_LocalAffineTransformVector.size(); t++)
    {
    this->m_TrajectoryWeightImageVector[t] = this->ComputeTrajectoryWeightImage(
      this->m_TrajectoryImageVector[t], this->m_BoundaryWeightImage);
    PicslImageHelper::WriteImage<WeightImageType>
      (this->m_TrajectoryWeightImageVector[t], "tmpTrajWeight.nii", t*10000+stopTime);
    }

  this->ComputeWeightedSumOfVelocityFields();
}

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
  return this->m_DisplacementField;
}

template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::ComputeDisplacementField()
{
  this->InitializeBuffers();
  
  typedef typename DisplacementFieldType::PixelType VectorType;
  VectorType disp;
  disp.Fill(0);
  this->m_DisplacementField->FillBuffer(disp); //identity transform

  int stopTime = -1;
  while (stopTime < this->m_TimeStampMax)
    {
    this->ComputeVelocityFieldBeforeOverlap(stopTime);
    PicslImageHelper::WriteDisplacementField<DisplacementFieldType>
      (this->m_VelocityField, "tmpVelocitySum.nii", stopTime);

    DisplacementFieldPointer exponentialField = 
      this->ComputeExponentialDisplacementField(this->m_VelocityField);
    PicslImageHelper::WriteDisplacementField<DisplacementFieldType>
      (exponentialField, "tmpExponentField.nii", stopTime);
 
    typedef itk::ComposeDisplacementFieldsImageFilter<DisplacementFieldType> ComposerType;
    typename ComposerType::Pointer composer = ComposerType::New();
    composer->SetDisplacementField( exponentialField );
    composer->SetWarpingField( this->m_DisplacementField );
    composer->Update();

    this->m_DisplacementField = composer->GetOutput();
    this->m_DisplacementField->DisconnectPipeline();
    }
}

template< class TScalarType,
          unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::DisplacementFieldPointer
PolyAffineTransform< TScalarType, NDimensions >
::ComputeExponentialDisplacementField(const DisplacementFieldPointer &velocityField)
{
  ExponentialImageFilterPointer filter = ExponentialImageFilterType::New();
  filter->SetInput(velocityField);

  filter->SetAutomaticNumberOfIterations(true);
  filter->Update();

  DisplacementFieldPointer exponentialField = filter->GetOutput();
  return exponentialField;
}

} // namespace

#endif
