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
#ifndef __itkLocalAffineTransform_h
#define __itkLocalAffineTransform_h

#include "itkImage.h"
#include "itkPointSet.h"
#include "itkSpatialObject.h"
#include "itkAffineTransform.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkPicslImageHelper.h"

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
class LocalAffineTransform:
  public AffineTransform< TScalarType, NDimensions >
{
public:
  /** Standard typedefs   */
  typedef LocalAffineTransform Self;
  typedef AffineTransform< TScalarType, NDimensions > Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(LocalAffineTransform, AffineTransform);

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

  typedef AffineTransform< TScalarType, NDimensions >   VelocityAffineTransformType;
  typedef typename VelocityAffineTransformType::Pointer          VelocityAffineTransformPointer;
  virtual VelocityAffineTransformType* GetVelocityAffineTransform();

  void AddFixedPoint(const InputPointType point);
  
  /** Set the Mask Image.  */
  itkSetObjectMacro(FixedMaskImage, MaskImageType);
  /** Get the Mask Image. */
  itkGetObjectMacro(FixedMaskImage, MaskImageType);
  /** Set the Mask Image.  */
  itkSetObjectMacro(MovingMaskImage, MaskImageType);
  /** Get the Mask Image. */
  itkGetObjectMacro(MovingMaskImage, MaskImageType);

  itkSetObjectMacro(DenseFixedPointSet, PointSetType);
  itkGetObjectMacro(DenseFixedPointSet, PointSetType);

  /** Set the Mask Image from a Spatial Object.  */
  template< class TSpatialObject > 
  void ComputeFixedMaskImageFromSpatialObject(TSpatialObject *spatialObject, SizeType size);
  template< class TSpatialObject > 
  void ComputeFixedMaskImageFromSpatialObject(TSpatialObject *spatialObject,
                                PointType origin,
                                DirectionType direction,
                                SpacingType spacing,
                                SizeType size);
  void ComputeMovingMaskImageAndDenseFixedPointSet();
  void AddMaskToPointSet(PointSetPointer &pointSet, const MaskImagePointer &mask);

protected:
  /** Construct an LocalAffineTransform object
   *
   * This method constructs a new LocalAffineTransform object and
   * initializes the matrix and offset parts of the transformation
   * to values specified by the caller.  If the arguments are
   * omitted, then the LocalAffineTransform is initialized to an identity
   * transformation in the appropriate number of dimensions.   */
  LocalAffineTransform(const MatrixType & matrix,
                  const OutputVectorType & offset);
  LocalAffineTransform(unsigned int paramDims);
  LocalAffineTransform();

  /** Destroy an LocalAffineTransform object   */
  virtual ~LocalAffineTransform();

  /** Print contents of an LocalAffineTransform */
  void PrintSelf(std::ostream & s, Indent indent) const;

  void ComputePrincipalLogorithm();

private:

  LocalAffineTransform(const Self & other);
  const Self & operator=(const Self &);

  VelocityAffineTransformPointer m_VelocityAffineTransform;
  PointSetPointer m_FixedPointSet;

  //points in the moving mask warped to fixed domain
  PointSetPointer m_DenseFixedPointSet;

  MaskImagePointer m_FixedMaskImage;
  MaskImagePointer m_MovingMaskImage;

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
#include "Templates/itkLocalAffineTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkLocalAffineTransform.hxx"
#endif

#endif /* __itkLocalAffineTransform_h */
