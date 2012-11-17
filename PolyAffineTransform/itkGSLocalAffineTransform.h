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
#ifndef __itkGSLocalAffineTransform_h
#define __itkGSLocalAffineTransform_h

#include "itkImage.h"
#include "itkPointSet.h"
#include "itkSpatialObject.h"
#include "itkAffineTransform.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkPicslImageHelper.h"
#include "itkTimeProbe.h"

namespace itk
{
/**
 * Affine transformation with local control points.
 *
 * This class serves PolyAffineTransform by attaching an impact region to
 * AffineTransform. It computes the principal logorithm of the homogenous
 * matrix which will be used for the velocity field. The impact region is
 * defined by a SpatialObject or an image mask.
 *
 * There are two template parameters for this class:
 *
 * ScalarT       The type to be used for scalar numeric values.  Either
 *               float or double.
 *
 * NDimensions   The number of dimensions of the vector space.
 *
 * \ingroup ITKTransform
 */

template<
  class TScalarType = double,      // Data type for scalars
                                   //    (e.g. float or double)
  unsigned int NDimensions = 3 >
// Number of dimensions in the input space
class GSLocalAffineTransform:
  public AffineTransform< TScalarType, NDimensions >
{
public:
  /** Standard typedefs   */
  typedef GSLocalAffineTransform Self;
  typedef AffineTransform< TScalarType, NDimensions > Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(GSLocalAffineTransform, AffineTransform);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Dimension of the domain space. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, NDimensions);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, NDimensions);
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
  itkStaticConstMacro( ParametersDimension, unsigned int,
                       NDimensions *( NDimensions + 1 ) );

  /** Parameters Type   */
  typedef typename Superclass::ParametersType            ParametersType;
  typedef typename Superclass::JacobianType              JacobianType;
  typedef typename Superclass::ScalarType                ScalarType;
  typedef typename Superclass::InputPointType            InputPointType;
  typedef typename Superclass::OutputPointType           OutputPointType;
  typedef typename Superclass::InputVectorType           InputVectorType;
  typedef typename Superclass::OutputVectorType          OutputVectorType;
  typedef typename Superclass::InputVnlVectorType        InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType       OutputVnlVectorType;
  typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;
  typedef typename Superclass::MatrixType                MatrixType;
  typedef typename Superclass::InverseMatrixType         InverseMatrixType;
  typedef typename Superclass::CenterType                CenterType;
  typedef typename Superclass::OffsetType                OffsetType;
  typedef typename Superclass::TranslationType           TranslationType; 

  /** Base inverse transform type. This type should not be changed to the
   * concrete inverse transform type or inheritance would be lost.*/
  typedef typename Superclass::InverseTransformBaseType InverseTransformBaseType;
  typedef typename InverseTransformBaseType::Pointer    InverseTransformBasePointer;

  typedef itk::PointSet<int, NDimensions>               PointSetType;
  typedef typename PointSetType::Pointer                PointSetPointer;

  typedef itk::Image<unsigned short, NDimensions>       MaskImageType;
  typedef typename MaskImageType::Pointer               MaskImagePointer;

  typedef typename MaskImageType::SizeType              SizeType;
  typedef typename MaskImageType::IndexType             IndexType;
  typedef typename MaskImageType::PointType             PointType;
  typedef typename MaskImageType::RegionType            RegionType;
  typedef typename MaskImageType::SpacingType           SpacingType;
  typedef typename MaskImageType::DirectionType         DirectionType;
  typedef typename itk::ContinuousIndex<TScalarType, NDimensions>
                                                        ContinuousIndexType;

  typedef Superclass                                    AffineTransformType;
  typedef typename AffineTransformType::Pointer         AffineTransformPointer;

  void AddFixedPoint(const InputPointType point);
  
  /** Set the Mask Image.  */
  itkSetObjectMacro(FixedMaskImage, MaskImageType);
  /** Get the Mask Image. */
  itkGetObjectMacro(FixedMaskImage, MaskImageType);
  /** Set the Mask Image.  */
  itkSetObjectMacro(MovingMaskImage, MaskImageType);
  /** Get the Mask Image. */
  itkGetObjectMacro(MovingMaskImage, MaskImageType);

  itkSetObjectMacro(SamplePointSet, PointSetType);
  itkGetObjectMacro(SamplePointSet, PointSetType);

  itkSetMacro(StartTime, int);
  itkGetMacro(StartTime, int);

  itkSetMacro(StopTime, int);
  itkGetMacro(StopTime, int);

  itkSetMacro(TimeStampMax, int);
  itkGetMacro(TimeStampMax, int);

  itkSetMacro(Overlapped, bool);
  itkGetMacro(Overlapped, bool);

  itkSetMacro(OverlapPointId, unsigned int);
  itkGetMacro(OverlapPointId, unsigned int);

  /** Methods to get timers */
  itkGetConstReferenceMacro(TimerComputeFixedMaskImage,  TimeProbe);
  itkGetConstReferenceMacro(TimerComputeMovingMaskImage, TimeProbe);
  itkGetConstReferenceMacro(TimerMatrixExponential,      TimeProbe);
  itkGetConstReferenceMacro(TimerMatrixLogarithm,        TimeProbe);
  itkGetConstReferenceMacro(TimerComputeSamplePointSet,  TimeProbe);
  itkGetConstReferenceMacro(TimerDilateFixedMaskImage,   TimeProbe);

  /** Compute the fixed mask image from a Spatial Object.  */
  template< class TSpatialObject > 
  void ComputeFixedMaskImageFromSpatialObject(TSpatialObject *spatialObject, SizeType size);

  /** Compute the fixed mask image from a Spatial Object with information 
   *  about origin, orientation and spacing.  */
  template< class TSpatialObject > 
  void ComputeFixedMaskImageFromSpatialObject(TSpatialObject *spatialObject,
                                PointType origin,
                                DirectionType direction,
                                SpacingType spacing,
                                SizeType size);

  /** Dilate the fixed mask image with a radius. */
  void DilateFixedMaskImage(unsigned int radius);

  /**
   * Transform the fixed mask into the moving mask. The moving mask
   * domain include both the fixed and moving domain.
   */
  void ComputeMovingMaskImage();

  /**
   * Compute the SamplePointSet at timeStamp. It warps the fixed mask
   * and the moving mask into a virtual domain at timeStamp. 
   */
  void ComputeSamplePointSet(int timeStamp);

  /** Compute the principal logarithm of the homgoeneous matrix of an affine transform. */
  vnl_matrix<TScalarType> ComputeLogarithmMatrix(AffineTransformPointer affineTransform);

  /** Compute the exponential transform from a homogeneous matrix. */
  AffineTransformPointer ComputeExponentialTransform(vnl_matrix<TScalarType> velocity);

  /**
   * The existing method AffineTransform::Scale(factor) scales around the center.
   * Instead, we scale around (0,0) in ScaleMatrixOffset().
   */
  AffineTransformPointer ScaleMatrixOffset(AffineTransformPointer velocity, double factor);

  /**
   * Convert a homogeneous matrix into an affine transform.
   */
  void SetHomogeneousTransform(AffineTransformPointer &affineTransform,
    const vnl_matrix<TScalarType> &homoMatrix);

  /** Get the homogeneous matrix from an affine transform. */
  void GetHomogeneousMatrix(vnl_matrix<TScalarType> &homoMatrix, 
    const AffineTransformPointer &affineTransform);

  /** Compute the partial transform during [startTime, stopTime].
   *  It is computed by exp(t*log(T)) where t=(stopTime-startTime)/m_TimeStampMax,
   *  and T is this transform. */
  AffineTransformPointer GetPartialTransform(int startTime, int stopTime);

  /** Compute the homogeneous matrix of the velocity of this transform. */
  vnl_matrix<TScalarType> GetVelocityMatrix();

  /** Get velocity field at the point **/
  OutputVectorType GetVelocityAtPoint(const PointType &point);


  /** Update the partial velocity matrix during current [m_StartTime, m_StopTime]. */
  void UpdatePartialVelocityMatrix();

  /** Get the partial velocity vector at a point during current [m_StartTime, m_StopTime]. */ 
  OutputVectorType GetPartialVelocityAtPoint(const PointType &point);
  
  /** Warp a mask into a PointSet by a transform. */
  void WarpMaskIntoPointSet(PointSetPointer &pointSet, const MaskImagePointer &mask,
                         const AffineTransformPointer &transform);
  /** Warp the fixed mask into a PointSet in a virtual domain at timeStamp. */
  void WarpFixedMaskIntoPointSet(int timeStamp);
  /** Warp the moving mask into a PointSet in a virtual domain at timeStamp. */
  void WarpMovingMaskIntoPointSet(int timeStamp);

protected:
  /** Construct an LocalAffineTransform object
   *
   * This method constructs a new LocalAffineTransform object and
   * initializes the matrix and offset parts of the transformation
   * to values specified by the caller.  If the arguments are
   * omitted, then the LocalAffineTransform is initialized to an identity
   * transformation in the appropriate number of dimensions.   */
  GSLocalAffineTransform(const MatrixType & matrix,
                  const OutputVectorType & offset);
  GSLocalAffineTransform(unsigned int paramDims);
  GSLocalAffineTransform();

  /** Destroy an LocalAffineTransform object   */
  virtual ~GSLocalAffineTransform();

  /** Print contents of an LocalAffineTransform */
  void PrintSelf(std::ostream & s, Indent indent) const;

private:

  GSLocalAffineTransform(const Self & other);
  const Self & operator=(const Self &);

  /** Each transform has its own starting time and stopping time of 
   *  its current trajectory. */
  int m_StartTime, m_StopTime;

  /** The number time stamps, usually it is 2^N with the least N such
   *  that 2^N >= max(size(image)), the maximum number of voxels on
   *  each dimension.
   */
  int m_TimeStampMax;

  bool m_Overlapped;
  unsigned int m_OverlapPointId;

  //The logarithm matrix of the homogeneous matrix of this transform .
  vnl_matrix<TScalarType> m_VelocityMatrix;

  //m_VelocityMatrix scaled by current time period.
  vnl_matrix<TScalarType> m_PartialVelocityMatrix;

  //sample points from which we compute trajectories.
  //now it contains points at fixed image domain.
  PointSetPointer m_SamplePointSet;


  MaskImagePointer m_FixedMaskImage;
  MaskImagePointer m_MovingMaskImage;

  itk::TimeProbe            m_TimerComputeFixedMaskImage;
  itk::TimeProbe            m_TimerComputeMovingMaskImage;
  itk::TimeProbe            m_TimerMatrixExponential;
  itk::TimeProbe            m_TimerMatrixLogarithm;
  itk::TimeProbe            m_TimerComputeSamplePointSet;
  itk::TimeProbe            m_TimerDilateFixedMaskImage;

}; //class LocalAffineTransform

}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_LocalAffineTransform(_, EXPORT, TypeX, TypeY)                \
  namespace itk                                                              \
  {                                                                          \
  _( 2 ( class EXPORT LocalAffineTransform< ITK_TEMPLATE_2 TypeX > ) )            \
  namespace Templates                                                        \
  {                                                                          \
  typedef LocalAffineTransform< ITK_TEMPLATE_2 TypeX >  LocalAffineTransform##TypeY; \
  }                                                                          \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkGSLocalAffineTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkGSLocalAffineTransform.hxx"
#endif

#endif /* __itkLocalAffineTransform_h */
