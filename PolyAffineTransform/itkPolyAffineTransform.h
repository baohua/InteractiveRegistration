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
#ifndef __itkPolyAffineTransform_h
#define __itkPolyAffineTransform_h

#include <iostream>

#include "itkDisplacementFieldTransform.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkExponentialDisplacementFieldImageFilter.h"
#include "itkComposeDisplacementFieldsImageFilter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkTimeProbe.h"

#include "itkLocalAffineTransform.h"
#include "itkPicslImageHelper.h"

namespace itk
{
/**
 * PolyAffineTransform is a diffeomorphic transform that combines several
 * local affine transforms. Each local affine transform is inherited from
 * itk::AffineTransform and is given an initial impact region designated
 * by a spatial object or an image mask.
 *
 * \ingroup Transforms
 *
 */

template< class TScalarType,
          unsigned int NDimensions = 3
>
// Number of dimensions in the output space
class PolyAffineTransform:
  public Transform< TScalarType, NDimensions, NDimensions >
{
public:
  /** Standard typedefs   */
  typedef PolyAffineTransform                   Self;
  typedef Transform< TScalarType,
                     NDimensions,
                     NDimensions >        Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(PolyAffineTransform, Transform);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Dimension of the domain space. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, NDimensions);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, NDimensions);
  //itkStaticConstMacro( ParametersDimension, unsigned int,
  //                     NDimensions * ( NDimensions + 1 ) * NTransforms);

  /** Parameters Type   */
  typedef typename Superclass::ParametersType      ParametersType;
  typedef typename Superclass::ParametersValueType ParametersValueType;

  /** Jacobian Type   */
  typedef typename Superclass::JacobianType JacobianType;

  /** Standard scalar type for this class */
  typedef typename Superclass::ScalarType ScalarType;

  /** Standard vector type for this class   */
  typedef Vector< TScalarType,
                  itkGetStaticConstMacro(InputSpaceDimension) >  InputVectorType;
  typedef Vector< TScalarType,
                  itkGetStaticConstMacro(OutputSpaceDimension) > OutputVectorType;
  typedef typename OutputVectorType::ValueType OutputVectorValueType;

  /** Standard covariant vector type for this class   */
  typedef CovariantVector< TScalarType,
                           itkGetStaticConstMacro(InputSpaceDimension) >
  InputCovariantVectorType;
  typedef CovariantVector< TScalarType,
                           itkGetStaticConstMacro(OutputSpaceDimension) >
  OutputCovariantVectorType;

  typedef typename Superclass::InputVectorPixelType  InputVectorPixelType;
  typedef typename Superclass::OutputVectorPixelType OutputVectorPixelType;

  /** Standard vnl_vector type for this class   */
  typedef vnl_vector_fixed< TScalarType,
                            itkGetStaticConstMacro(InputSpaceDimension) >
  InputVnlVectorType;
  typedef vnl_vector_fixed< TScalarType,
                            itkGetStaticConstMacro(OutputSpaceDimension) >
  OutputVnlVectorType;

  /** Standard coordinate point type for this class   */
  typedef Point< TScalarType,
                 itkGetStaticConstMacro(InputSpaceDimension) >
  InputPointType;
  typedef typename InputPointType::ValueType InputPointValueType;
  typedef Point< TScalarType,
                 itkGetStaticConstMacro(OutputSpaceDimension) >
  OutputPointType;
  typedef typename OutputPointType::ValueType OutputPointValueType;

  /** Standard matrix type for this class   */
  typedef Matrix< TScalarType, itkGetStaticConstMacro(OutputSpaceDimension),
                  itkGetStaticConstMacro(InputSpaceDimension) >
  MatrixType;
  typedef typename MatrixType::ValueType MatrixValueType;

  /** Standard inverse matrix type for this class   */
  typedef Matrix< TScalarType, itkGetStaticConstMacro(InputSpaceDimension),
                  itkGetStaticConstMacro(OutputSpaceDimension) >
  InverseMatrixType;

  typedef InputPointType CenterType;

  typedef OutputVectorType               OffsetType;
  typedef typename OffsetType::ValueType OffsetValueType;

  typedef OutputVectorType TranslationType;

  typedef typename TranslationType::ValueType TranslationValueType;

  /** Base inverse transform type. This type should not be changed to the
   * concrete inverse transform type or inheritance would be lost. */
  typedef typename Superclass::InverseTransformBaseType InverseTransformBaseType;
  typedef typename InverseTransformBaseType::Pointer    InverseTransformBasePointer;

  /** Types for local transforms */
  typedef LocalAffineTransform<TScalarType, NDimensions>       LocalAffineTransformType;
  typedef typename LocalAffineTransformType::Pointer           LocalAffineTransformPointer;
  typedef typename LocalAffineTransformType::ParametersType    LocalAffineParametersType;
  typedef typename LocalAffineTransformType::PointSetType      PointSetType;
  typedef typename LocalAffineTransformType::PointSetPointer   PointSetPointer;
  typedef typename LocalAffineTransformType::AffineTransformType 
                                                               AffineTransformType;
  typedef typename LocalAffineTransformType::AffineTransformPointer
                                                               AffineTransformPointer;

  typedef itk::DisplacementFieldTransform<TScalarType, NDimensions> 
                                                               DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::Pointer     DisplacementFieldTransformPointer;

  typedef typename DisplacementFieldTransformType::DisplacementFieldType
                                                               DisplacementFieldType;
  typedef typename DisplacementFieldType::Pointer              DisplacementFieldPointer;    
  typedef typename DisplacementFieldType::PixelType            DisplacementVectorType;

  typedef typename LocalAffineTransformType::MaskImageType     MaskImageType;
  typedef typename MaskImageType::Pointer                      MaskImagePointer;
  
  typedef typename MaskImageType::IndexType                    IndexType;
  typedef typename itk::ContinuousIndex<TScalarType, NDimensions>
                                                               ContinuousIndexType;
  typedef typename MaskImageType::PointType                    PointType;
  typedef typename MaskImageType::SizeType                     SizeType;
  typedef typename MaskImageType::RegionType                   RegionType;
  typedef typename MaskImageType::SpacingType                  SpacingType;
  typedef typename MaskImageType::DirectionType                DirectionType;

  typedef typename itk::Image<TScalarType, NDimensions>        DistanceMapImageType;
  typedef typename DistanceMapImageType::Pointer               DistanceMapImagePointer;

  typedef std::vector<PointSetPointer>                         PointSetVectorType;

  /** TrajectoryImageType: pixel values denotes time stamps */
  typedef itk::Image<int, NDimensions>                         TrajectoryImageType;
  typedef typename TrajectoryImageType::Pointer                TrajectoryImagePointer;

  typedef itk::ExponentialDisplacementFieldImageFilter
    <DisplacementFieldType, DisplacementFieldType>             ExponentialImageFilterType;
  typedef typename ExponentialImageFilterType::Pointer         ExponentialImageFilterPointer;

  /** Transform queue type */
  typedef std::vector<LocalAffineTransformPointer>             LocalAffineTransformVectorType;
  typedef std::vector<TrajectoryImagePointer>                  TrajectoryImageVectorType;
  typedef std::vector<DistanceMapImagePointer>                 DistanceMapImageVectorType;

  LocalAffineTransformVectorType & GetLocalAffineTransformVector()
    {
    return m_LocalAffineTransformVector;
    }
  unsigned int GetNumberOfLocalAffineTransforms()
    {
    return m_LocalAffineTransformVector.size();
    }

  /** Set the transformation to an Identity
   *
   * This sets the matrix to identity and the Offset to null. */
  virtual void SetIdentity(void);

  /** Set the transformation from a container of parameters.
   * The first (NOutputDimension x NInputDimension) parameters define the
   * matrix and the last NOutputDimension parameters the translation.
   * Offset is updated based on current center. */
  void SetParameters(const ParametersType & parameters);

  /** Get the Transformation Parameters. */
  const ParametersType & GetParameters(void) const;

  /** Set the fixed parameters and update internal transformation. */
  virtual void SetFixedParameters(const ParametersType &);

  /** Get the Fixed Parameters. */
  virtual const ParametersType & GetFixedParameters(void) const;

  /** Transform by an affine transformation
   *
   * This method applies the affine transform given by self to a
   * given point or vector, returning the transformed point or
   * vector.  The TransformPoint method transforms its argument as
   * an affine point, whereas the TransformVector method transforms
   * its argument as a vector. */
  OutputPointType       TransformPoint(const InputPointType & point) const;

  OutputVectorType      TransformVector(const InputVectorType & vector) const;

  OutputVectorType      TransformVector(const InputVectorType & vector,
                                        const InputPointType & itkNotUsed(point) ) const
    { return TransformVector( vector ); }

  OutputVnlVectorType   TransformVector(const InputVnlVectorType & vector) const;

  OutputVnlVectorType   TransformVector(const InputVnlVectorType & vector,
                                        const InputPointType & itkNotUsed(point) ) const
    { return TransformVector( vector ); }

  //OutputVectorPixelType TransformVector(const InputVectorPixelType & vector) const;

  //OutputVectorPixelType TransformVector(const InputVectorPixelType & vector,
  //                                      const InputPointType & itkNotUsed(point) ) const
  //  { return TransformVector( vector ); }


  OutputCovariantVectorType TransformCovariantVector(
      const InputCovariantVectorType & vector) const;

  OutputCovariantVectorType TransformCovariantVector(
      const InputCovariantVectorType & vector,
      const InputPointType & itkNotUsed(point) ) const
    { return TransformCovariantVector( vector ); }

  //OutputVectorPixelType TransformCovariantVector(
  //    const InputVectorPixelType & vector) const;

  //OutputVectorPixelType TransformCovariantVector(
  //    const InputVectorPixelType & vector,
  //    const InputPointType & itkNotUsed(point) ) const
  //  { return TransformCovariantVector( vector ); }

  /** Compute the Jacobian of the transformation
   *
   * This method computes the Jacobian matrix of the transformation.
   * given point or vector, returning the transformed point or
   * vector. The rank of the Jacobian will also indicate if the transform
   * is invertible at this point. */
  const JacobianType & GetJacobian(const InputPointType & point) const;

  /** get local Jacobian for the given point
   * \c j will sized properly as needed.
   * This is a thread-safe version for GetJacobian(). Otherwise,
   * m_Jacobian could be changed for different values in different threads. */
  void GetJacobianWithRespectToParameters(const InputPointType  &x,
                                          JacobianType &j) const;

  /** Create inverse of an affine transformation
   *
   * This populates the parameters an affine transform such that
   * the transform is the inverse of self. If self is not invertible,
   * an exception is thrown.
   * Note that by default the inverese transform is centered at
   * the origin. If you need to compute the inverse centered at a point, p,
   *
   * \code
   * transform2->SetCenter( p );
   * transform1->GetInverse( transform2 );
   * \endcode
   *
   * transform2 will now contain the inverse of transform1 and will
   * with its center set to p. Flipping the two statements will produce an
   * incorrect transform.
   *
   */
  bool GetInverse(Self *inverse) const;

  /** Return an inverse of this transform. */
  virtual InverseTransformBasePointer GetInverseTransform() const;

  itkSetMacro(TimeStampAutomaticFlag, bool);
  itkGetMacro(TimeStampAutomaticFlag, bool);

  itkSetMacro(PadTrajectory, int);
  itkGetMacro(PadTrajectory, int);
  
  itkSetMacro(PadBoundary, int);
  itkGetMacro(PadBoundary, int);

  itkSetMacro(StopAllTrajectoriesAtOverlap, bool);
  itkGetMacro(StopAllTrajectoriesAtOverlap, bool);

  itkSetMacro(DecayRateOfBoundary, double);
  itkGetMacro(DecayRateOfBoundary, double);

  itkSetMacro(DecayRateOfTrajectory, double);
  itkGetMacro(DecayRateOfTrajectory, double);

  void AddLocalAffineTransform( LocalAffineTransformType *t  )
  {
    LocalAffineTransformPointer tp = t;
    m_LocalAffineTransformVector.push_back(tp);
    this->Modified();
  }

  void SetTimeStampLog(unsigned int timeStampLog)
  {
    this->m_TimeStampLog = timeStampLog;

    //m_TimeStampMax is the number of time stamps
    this->m_TimeStampMax = 1 << timeStampLog;
  }

  itkGetMacro(TimeStampLog, int);
  itkGetMacro(TimeStampMax, int);

  void UpdateTimeStampIfAutomatic();

  /** Translate the time stamp to a nonegative number to be stored
   *  in a trajectory image. The trajectory background uses 0.
   */
  int TimeStampToImageValue(int timeStep)
  { 
    //background is 0
    return timeStep + 1;
  }
  int ImageValueToTimeStamp(int imageValue)
  { 
    return imageValue - 1;
  }

  /** Get the distance between voxel at index and at index+1.
   *  index+1 is for the next voxel with all index[i] added by 1.
   */
  TScalarType GetDiagonalSpacing(MaskImagePointer mask);

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
  bool InitializeBuffers();

  /** Allocate memory for DisplacementField and fill it with zeros.
   */
  void InitializeDisplacementField();

  /** Compute the boundary as the union of all moving and fixed image
   *  domains of local transforms.
   */
  void InitializeBoundaryMask();

  /** Compute the minimum stop time stamps of all local transforms. 
   *  It also prints the stop time stamps if debugInfo is not null.
   */
  int GetMinStopTime(char *debugInfo = NULL);

  /** Rewind a trajectory of a transform back to stopTime.
   */
  void RewindTrajectory(unsigned int transformId, int stopTime);

  /** Move the current trajectory of a transform one step forward.
   */
  void ComputeNextStepTrajectory(unsigned int transformId);

  /** Check if the point specified by index exists in trajectories of 
   *  transforms other than the transform specified by transformId.
   */
  bool PointExistsInOtherTrajectories(unsigned int transformId, IndexType index);

  /** Combine multiple trajectories into one trajectory. 
   */
  void CombineTrajectories();

  /** Compute the distance map of the boundary mask.
   */
  DistanceMapImagePointer ComputeBoundaryDistanceMapImage();

  /** Compute the distance map of a trajectory.
   */
  DistanceMapImagePointer ComputeTrajectoryDistanceMapImage(TrajectoryImagePointer traj);

  /** Comopute the weighted sum of velocity fields specified by the local
   *  affine transforms.
   */
  void ComputeWeightedSumOfVelocityFields();

  /** Compute the weighted sum of velocity vectors at a point with given
   *  distances to different trajectories and the boundary. 
   */
  void ComputeWeightedSumOfVelocitiesAtPoint(DisplacementVectorType &velocitySum, const PointType &point,
    const double distances[], double distanceToCombinedTraj, double distanceToImageBoundary);

  /** Compute the exponential mapping a velocity field.
   */
  DisplacementFieldPointer ComputeExponentialDisplacementField(
    const DisplacementFieldPointer &velocityField);

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
  void ComputeVelocityFieldBeforeOverlap();

  /** Compute the displacement field of this PolyAffineTransform.
   *  (1) It initializes the trajectories of local transforms.
   *  (2) It will compute the velocity sum for each segment of time
   *  during which there is no trajectory overlap. 
   *  (3) It computes the exponential mapping to get a displacement field
   *  during that time segment.
   *  (4) It computes the composition of these displacement fields at
   *  mutiple tiem segments.
   */
  void ComputeDisplacementField();

  /** Get the displacement field of this PolyAffineTransform.
   */
  virtual DisplacementFieldType* GetDisplacementField();

  /** Print the timers. */   
  void PrintTimers();

protected:
  /** Construct an PolyAffineTransform object
   *
   * This method constructs a new PolyAffineTransform object and
   * initializes the matrix and offset parts of the transformation
   * to values specified by the caller.  If the arguments are
   * omitted, then the PolyAffineTransform is initialized to an identity
   * transformation in the appropriate number of dimensions. */
  PolyAffineTransform();

  /** Destroy an PolyAffineTransform object */
  virtual ~PolyAffineTransform();

  /** Print contents of an PolyAffineTransform */
  void PrintSelf(std::ostream & s, Indent indent) const;

private:

  PolyAffineTransform(const Self & other);
  const Self & operator=(const Self &);

  // Flag to set m_TimeStampMax and m_TimeStampLog automatically.
  // If it is true, it will overwrite the user's given values.
  bool                                      m_TimeStampAutomaticFlag;

  // m_TimeStampLog is the logarithm of the number of time stamps
  int                                       m_TimeStampLog;

  /** The number time stamps, usually it is 2^N with the least N such
   *  that 2^N >= max(size(image)), the maximum number of voxels on
   *  each dimension.
   */
  int                                       m_TimeStampMax;

  // Pad to the image boundary
  int                                       m_PadBoundary;
  // Radius to dilate the trajectory
  int                                       m_PadTrajectory;
  // A switch to decide if it stops all trajectories at overlap
  bool                                      m_StopAllTrajectoriesAtOverlap;

  // For exponential decay rate of boundary distance
  double                                    m_DecayRateOfBoundary;
  // For exponential decay rate of trajectory distance
  double                                    m_DecayRateOfTrajectory;

  // Array of the local affine transforms
  LocalAffineTransformVectorType            m_LocalAffineTransformVector;  

  // The velocity field summed from different local transforms
  DisplacementFieldPointer                  m_VelocityField;
  // The total displacement field
  DisplacementFieldPointer                  m_DisplacementField;

  // Boundary mask image includes the all image domains of local transforms
  // It also apply a padding around it border specified by m_PadBoundary.
  MaskImagePointer                          m_BoundaryMask;

  /** Array of the trajectories of local transforms. m_TrajectoryImageVector[t]
   *  stores the t-th trajectory during [m_StartTime, m_StopTime] of the t-th
   *  local transform.
   */
  TrajectoryImageVectorType                 m_TrajectoryImageVector;  

  // Combined image of the trajectories of local transforms 
  TrajectoryImagePointer                    m_CombinedTrajectoryImage;

  // Distance map of the boundary image
  DistanceMapImagePointer                   m_BoundaryDistanceMapImage;
  // Array of distance maps of the trajectories
  DistanceMapImageVectorType                m_TrajectoryDistanceMapImageVector;
  // Distance map to the combined trajectory
  DistanceMapImagePointer                   m_CombinedTrajectoryDistanceMapImage;

  // Timers for different methods
  itk::TimeProbe                            m_TimerComputeNextStepTrajectory;
  itk::TimeProbe                            m_TimerComputeTrajectoryDistanceMapImage;
  itk::TimeProbe                            m_TimerComputeBoundaryDistanceMapImage;

  itk::TimeProbe                            m_TimerRewindTrajectory;
  itk::TimeProbe                            m_TimerCombineTrajectories;
  itk::TimeProbe                            m_TimerComputeWeightedSumOfVelocityFields;

  itk::TimeProbe                            m_TimerInitializeBuffers;
  itk::TimeProbe                            m_TimerInitializeBoundaryMask;
  itk::TimeProbe                            m_TimerInitializeDisplacementField;

  itk::TimeProbe                            m_TimerInitializeIteration;
  itk::TimeProbe                            m_TimerComputeVelocityFieldBeforeOverlap;
  itk::TimeProbe                            m_TimerExponentialMapping;
  itk::TimeProbe                            m_TimerDisplacementFieldComposing;

  itk::TimeProbe                            m_TimerComputeDisplacementField;

}; //class PolyAffineTransform
}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_PolyAffineTransform(_, EXPORT, TypeX, TypeY)                         \
  namespace itk                                                                                 \
  {                                                                                             \
  _( 3 ( class EXPORT PolyAffineTransform< ITK_TEMPLATE_3 TypeX > ) )                     \
  namespace Templates                                                                           \
  {                                                                                             \
  typedef PolyAffineTransform< ITK_TEMPLATE_3 TypeX > PolyAffineTransform##TypeY; \
  }                                                                                             \
  }

#if ITK_TEMPLATE_EXPLICIT
//template < class TScalarType, unsigned int NDimensions, unsigned int
// NDimensions>
//   const unsigned int itk::PolyAffineTransform< TScalarType,
// NDimensions, NDimensions >::InputSpaceDimension;
//template < class TScalarType, unsigned int NDimensions, unsigned int
// NDimensions>
//   const unsigned int itk::PolyAffineTransform< TScalarType,
// NDimensions, NDimensions >::OutputSpaceDimension;
//template < class TScalarType, unsigned int NDimensions, unsigned int
// NDimensions>
//   const unsigned int itk::PolyAffineTransform< TScalarType,
// NDimensions, NDimensions >::ParametersDimension;
#include "Templates/itkPolyAffineTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkPolyAffineTransform.hxx"
#endif

#endif /* __itkPolyAffineTransform_h */
