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
  this->m_DecayConstant = 1.0;
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
 * Initialize the trajectory of each transform at time stamp 0 according to
 * the inital spatial object for the local affine transform.
 */
template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::AddMaskToTrajectory(TrajectoryImagePointer &traj, 
                      const MaskImagePointer &mask, 
                      typename TrajectoryImageType::PixelType label)
{
  IndexType index, index2;
  PointType point;
  bool inside;
  
  typedef ImageRegionIteratorWithIndex< TrajectoryImageType > IteratorType;
  IteratorType it( traj, traj->GetLargestPossibleRegion() );

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    traj->TransformIndexToPhysicalPoint( index, point );
    inside = mask->TransformPhysicalPointToIndex( point, index2 );

    if( inside )
      {
      if (mask->GetPixel( index2 ) > 0)
        {
        it.Set(label); //forebround
        }
      }
    } //for iterator
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
::InitializeFrontier(FrontierType &frontier, const PointSetPointer &pointSet)
{
  SizeValueType num = pointSet->GetNumberOfPoints();
  frontier.resize(num);

  typename PointSetType::PointsContainer *container = pointSet->GetPoints();
  typename PointSetType::PointsContainer::ConstIterator it = container->Begin();
  while (it != container->End())
    {
    frontier.push_back(it->Value());
    PointType point = it->Value();
    if (point[0] <= 4 && point[1] <= 4) {
      int len = frontier.size();
      std::cout << "InitializeFrontier " << point << " y=" << frontier[len-1] << std::endl;
      }
    it++;
    }
}

/**
 * Compute the trajectory of the spatial object under each transform by
 * propagation at each time stamp.
 */
template< class TScalarType,
          unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::TrajectoryImagePointer
PolyAffineTransform< TScalarType, NDimensions >
::ComputeTrajectory(unsigned int transformId)
{
  int N = this->m_TimeStampNumber;

  LocalAffineTransformPointer trans = this->m_LocalAffineTransformVector[transformId];

  typedef typename LocalAffineTransformType::VelocityAffineTransformPointer VelocityAffineTransformPointer;
  VelocityAffineTransformPointer velocity = trans->GetVelocityAffineTransform();

  TrajectoryImagePointer traj = TrajectoryImageType::New();
  traj->CopyInformation(this->m_BoundaryMask);
  traj->SetRegions(this->m_BoundaryMask->GetLargestPossibleRegion());
  traj->Allocate();
  traj->FillBuffer(0);

  //this->AddMaskToTrajectory(traj, trans->GetFixedMaskImage(), 1);
  //this->AddMaskToTrajectory(traj, trans->GetMovingMaskImage(), 1);
  
  typedef itk::NearestNeighborInterpolateImageFunction< TrajectoryImageType, TScalarType >
                                                       InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(traj);

  //queueNext stores the unprocessed points in
  //trajectory(j/N, M) - trajectory((j-1)/N, M).
	FrontierType frontier0, frontier1;
  FrontierType *frontierNow, *frontierNext;
  this->InitializeFrontier(frontier0, trans->GetDenseFixedPointSet());

  MaskImagePointer domain = this->m_BoundaryMask;
	for (int ts=0; ts<=N; ts++)
    {
    if (ts % 2 == 0) 
      {
		  frontierNow = &frontier0;
		  frontierNext = &frontier1;
      }
    else
      {
		  frontierNow = &frontier1;
		  frontierNext = &frontier0;
      }
    frontierNext->clear();
    //std::cout << "frontNow size = " << frontierNow->size() << std::endl;

		for (int f=0; f<frontierNow->size(); f++)
      {
			PointType y = frontierNow->at(f);
      PointType z, step;
			//It will avoid accumulated error if we compute z=(1+V/N)*y by;
      //x=T^-1 (ts,y);
      //z=T (ts+1,x);
      step = velocity->TransformPoint( y );
      for (unsigned int d=0; d<NDimensions; d++)
        {
        z[d] = y[d] + step[d] / N;
        }
      frontierNext->push_back(z);
      //this->AddSegmentIntoTrajectory(y, z, traj, frontierNext, 0, ts+1);

      IndexType index;
      bool insideImage = traj->TransformPhysicalPointToIndex(y, index);
      if (insideImage)
        {
        if (traj->GetPixel(index) <= 0)
          {
          traj->SetPixel(index, ts+1); //timestamp + 1
          if (ts==0 && index[0] == 0 && index[1] == 0)
            std::cout << "add traj: y = " << y << std::endl;
          }
        }
      }
    } //for timestamp 0..N

  frontier0.clear();
  frontier1.clear();

  char fname[256];
  sprintf(fname, "tmpTraj%d.nii", transformId);
  itk::PicslImageHelper::WriteImage<TrajectoryImageType>(traj, fname);

  return traj;
}

/**
 * Given a startPoint and its next move to endPoint. Put all points on the
 * segment (startPoint, endPoint] to the trajectory and frontier.
 */
template< class TScalarType,
          unsigned int NDimensions >
void
PolyAffineTransform< TScalarType, NDimensions >
::AddSegmentIntoTrajectory(PointType startPoint, PointType endPoint,
                           TrajectoryImagePointer traj, FrontierType *frontier,
                           int background, int foreground)
{
  ContinuousIndexType middleIndex, startIndex, endIndex;
  typename ContinuousIndexType::VectorType diffIndex;
  PointType middlePoint;
  IndexType index;
  RegionType region;
  bool insideImage;
  TScalarType maxDiff;
  unsigned int maxDim;

  traj->TransformPhysicalPointToContinuousIndex(startPoint, startIndex);
  traj->TransformPhysicalPointToContinuousIndex(endPoint, endIndex);

  diffIndex = endIndex - startIndex;
  maxDim = 0;
  maxDiff = vcl_abs(diffIndex[0]);
  for (unsigned int d=1; d<NDimensions; d++)
    {
      if (maxDiff < vcl_abs(diffIndex[d]))
        {
        maxDim = d;
        maxDiff = vcl_abs(diffIndex[d]);
        }
    }

  // We need to interpolate the segment along the maxDim-th dimension.
  typename ContinuousIndexType::VectorType stepIndex = diffIndex / maxDiff;
  region = traj->GetLargestPossibleRegion();

  middleIndex = endIndex;

  for (int s=0; s<maxDiff; s++)
    {
    index.CopyWithRound(middleIndex);
    insideImage = region.IsInside(index);
    if (insideImage && traj->GetPixel(index) == background)
      {
      traj->SetPixel(index, foreground); //foreground = timestamp + 1
      traj->TransformContinuousIndexToPhysicalPoint(middleIndex, middlePoint); 
      frontier->push_back(middlePoint);      
      }
    for (unsigned int d=0; d<NDimensions; d++)
      {
      middleIndex[d] = middleIndex[d] - stepIndex[d];      
      }
    break;
    }
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
void
PolyAffineTransform< TScalarType, NDimensions >
::ComputeVelocityField()
{
  if (m_LocalAffineTransformVector.size() == 0) 
    {
    return;
    }

  this->InitializeBoundaryMask();
  PicslImageHelper::WriteImage<MaskImageType>
    (this->m_BoundaryMask, "tmpBoundMask.nii");

  unsigned int transformNumber = this->m_LocalAffineTransformVector.size();
  this->m_TrajectoryWeightImageVector.resize(transformNumber);
  this->m_TrajectoryImageVector.resize(transformNumber);

  this->m_BoundaryWeightImage = 
    this->ComputeBoundaryWeightImage();
  PicslImageHelper::WriteImage<WeightImageType>
    (this->m_BoundaryWeightImage, "tmpBoundWeight.nii");

  for (unsigned int t=0; t<transformNumber; t++)
    {
    this->m_TrajectoryImageVector[t] = this->ComputeTrajectory(t);
    this->m_TrajectoryWeightImageVector[t] =
      this->ComputeTrajectoryWeightImage(this->m_TrajectoryImageVector[t],
                                         this->m_BoundaryWeightImage);
    char fname[256];
    sprintf(fname, "tmpTrajWeight%d.nii", t);
    PicslImageHelper::WriteImage<WeightImageType>
      (this->m_TrajectoryWeightImageVector[t], fname);
    }

  this->m_VelocityField = DisplacementFieldType::New();
  this->m_VelocityField->CopyInformation(this->m_BoundaryMask);
  this->m_VelocityField->SetRegions(this->m_BoundaryMask->GetLargestPossibleRegion());
  this->m_VelocityField->Allocate();

  typedef typename LocalAffineTransformType::VelocityAffineTransformType 
    VelocityAffineTransformType;
  typedef typename LocalAffineTransformType::VelocityAffineTransformPointer
    VelocityAffineTransformPointer;

  PointType point, velocity;
  typename DisplacementFieldType::PixelType velocitySum;
  IndexType index;
  TScalarType weightTraj;

  typedef ImageRegionIteratorWithIndex< DisplacementFieldType > IteratorType;
  IteratorType it( this->m_VelocityField, this->m_VelocityField->GetLargestPossibleRegion() );

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    this->m_VelocityField->TransformIndexToPhysicalPoint( index, point );

    velocitySum.Fill(0.0);
    for (int i=0; i<transformNumber; i++)
      {
      LocalAffineTransformPointer localAffine = this->m_LocalAffineTransformVector[i];
      VelocityAffineTransformPointer localVelocity = localAffine->GetVelocityAffineTransform();

      velocity = localVelocity->TransformPoint(point);
      weightTraj = this->m_TrajectoryWeightImageVector[i]->GetPixel(index);

      for (unsigned int j=0; j<NDimensions; j++)
        {
        velocitySum[j] += weightTraj * velocity[j];
        }
      }
    it.Set(velocitySum);
    }
  PicslImageHelper::WriteDisplacementField<DisplacementFieldType>
    (this->m_VelocityField, "tmpVelocitySum.nii");

}

template< class TScalarType,
          unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::DisplacementFieldType*
PolyAffineTransform< TScalarType, NDimensions >
::GetVelocityField()
{
  if (this->m_VelocityField.IsNull())
    {
    this->ComputeVelocityField();
    }
  return this->m_VelocityField.GetPointer();
}

template< class TScalarType,
          unsigned int NDimensions >
typename PolyAffineTransform< TScalarType, NDimensions >::DisplacementFieldType*
PolyAffineTransform< TScalarType, NDimensions >
::GetExponentialDisplacementField()
{
  if (this->m_DisplacementField.IsNull())
    {
    ExponentialImageFilterPointer filter = ExponentialImageFilterType::New();
    filter->SetInput(this->GetVelocityField());

    filter->SetAutomaticNumberOfIterations(true);
    //filter->SetMaximumNumberOfIterations( this->m_TimeStampLog );
    filter->Update();

    this->m_DisplacementField = filter->GetOutput();
    }

  return this->m_DisplacementField.GetPointer();
}

} // namespace

#endif
