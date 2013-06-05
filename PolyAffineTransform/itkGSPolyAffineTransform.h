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


/*
 * itkGSPolyAffineTransform.h
 *
 *  Created on: Nov 14, 2012
 *      Author: songgang
 */

#ifndef __itkGSPolyAffineTransform_h
#define __itkGSPolyAffineTransform_h


#include "itkMultiTransformBase.h"
#include "itkConstantVelocityFieldTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkGSLocalAffineTransform.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkDanielssonDistanceMapImageFilter.h"

#include "itkImageFileWriter.h"

#include <vector>
#include "vnl/vnl_math.h"

#include "itkCompositeTransform.h"

namespace itk
{

/** \class GSPolyAffineTransform
 * \brief PolyAffine Transform Parameterization:
 *  Key components:
 *  1. Trajectory Collision:   a collision time list
 *  2. vectors for each collision period
 *    velocity field
 *    affine transform
 *    affine velocity transform
 *    update transform
 *
 *
 *
 *
 * Provides local/dense/high-dimensionality transformation via a
 * a velocity field.
 *
 * \author Gang Song
 *
 * \ingroup ITKDisplacementField
 */
template
<class TScalar, unsigned int NDimensions = 3>
class ITK_EXPORT GSPolyAffineTransform :
  public MultiTransformBase<TScalar, NDimensions, NDimensions>
{
public:
  /** Standard class typedefs. */
  typedef GSPolyAffineTransform                            Self;
  typedef MultiTransformBase<TScalar, NDimensions, NDimensions> Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( GSPolyAffineTransform, Transform );

  /** New macro for creation of through a Smart Pointer */
  itkNewMacro( Self );

  /** Component transform type **/
  typedef Superclass    TransformType;
  typedef typename Superclass::Pointer TransformTypePointer;
  /** InverseTransform type. */
  typedef typename Superclass::InverseTransformBasePointer  InverseTransformBasePointer;


  /** generic types for transform type **/

  /** Scalar type. */
  typedef typename Superclass::ScalarType          ScalarType;
  /** Parameters type. */
  typedef typename Superclass::ParametersType      ParametersType;
  typedef typename Superclass::ParametersValueType ParametersValueType;
  /** Derivative type */
  typedef typename Superclass::DerivativeType DerivativeType;
  /** Jacobian type. */
  typedef typename Superclass::JacobianType JacobianType;
  /** Transform category type. */
  typedef typename Superclass::TransformCategoryType TransformCategoryType;

  /* Types relative to the container transform. */

  /** Standard coordinate point type for this class. */
  typedef typename Superclass::InputPointType             InputPointType;
  typedef typename Superclass::OutputPointType            OutputPointType;
  /** Standard vector type for this class. */
  typedef typename Superclass::InputVectorType            InputVectorType;
  typedef typename Superclass::OutputVectorType           OutputVectorType;
  /** Standard covariant vector type for this class */
  typedef typename Superclass::InputCovariantVectorType   InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType  OutputCovariantVectorType;
  /** Standard vnl_vector type for this class. */
  typedef typename Superclass::InputVnlVectorType         InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType        OutputVnlVectorType;
  /** Standard Vectorpixel type for this class */
  typedef typename Superclass::InputVectorPixelType       InputVectorPixelType;
  typedef typename Superclass::OutputVectorPixelType      OutputVectorPixelType;
  /** Standard DiffusionTensor3D typedef for this class */
  typedef typename Superclass::InputDiffusionTensor3DType  InputDiffusionTensor3DType;
  typedef typename Superclass::OutputDiffusionTensor3DType OutputDiffusionTensor3DType;
  /** Standard SymmetricSecondRankTensor typedef for this class */
  typedef typename Superclass::InputSymmetricSecondRankTensorType InputSymmetricSecondRankTensorType;
  typedef typename Superclass::OutputSymmetricSecondRankTensorType  OutputSymmetricSecondRankTensorType;


  /** Transform queue type */
  typedef typename Superclass::TransformQueueType         TransformQueueType;


  /** The number of parameters defininig this transform. */
  typedef typename Superclass::NumberOfParametersType NumberOfParametersType;

  /** Dimension of the domain spaces. */
  itkStaticConstMacro( InputDimension, unsigned int, NDimensions );
  itkStaticConstMacro( OutputDimension, unsigned int, NDimensions );

  /** If all sub-transforms are linear, then the multi-transform is linear. */
  // virtual bool IsLinear() const;

  /** If all sub-transforms are of the same category, return that category.
   * Otherwise return UnknownTransformCategory. */
  // virtual TransformCategoryType GetTransformCategory() const;

  // /** Get/Set Parameter functions work on all sub-transforms.
  //     The parameter data from each sub-transform is
  //     concatenated into a single ParametersType object.
  //     \note The sub-transforms are read in forward queue order,
  //     so the returned array is ordered in the same way. That is,
  //     first sub-transform to be added is returned first in the
  //     parameter array.*/
  using Superclass::GetParameters;
  // virtual const ParametersType & GetParameters(void) const;

  // /* SetParameters for all sub-transforms.
  //  * See GetParameters() for parameter ordering. */
  using Superclass::SetParameters;
  // virtual void  SetParameters(const ParametersType & p);

  // /* GetFixedParameters for all sub-transforms.
  //  * See GetParameters() for parameter ordering. */
  using Superclass::GetFixedParameters;
  // virtual const ParametersType & GetFixedParameters(void) const;

  //  SetFixedParameters for all sub-transforms.
  //  * See GetParameters() for parameter ordering.
  using Superclass::SetFixedParameters; 
  // virtual void SetFixedParameters(const ParametersType & fixedParameters);

  // /* Get total number of parameters. Sum of all sub-transforms. */
  using Superclass::GetNumberOfParameters;
  // virtual NumberOfParametersType GetNumberOfParameters(void) const;

  // /* Get total number of local parameters, the sum of all sub-transforms. */
  using Superclass::GetNumberOfLocalParameters;
  // virtual NumberOfParametersType GetNumberOfLocalParameters(void) const;

  // /* Get total number of fixed parameters, the sum of all sub-transforms. */
  using Superclass::GetNumberOfFixedParameters;
  // virtual NumberOfParametersType GetNumberOfFixedParameters(void) const;

  // /** Update the transform's parameters by the values in \c update.
  //  * See GetParameters() for parameter ordering. */
  // virtual void UpdateTransformParameters( const DerivativeType & update, ScalarType  factor = 1.0 );

  // virtual bool GetInverse( Self *inverse ) const;



  /** Compute the position of point in the new space.
  *
  * Transforms are applied starting from the *back* of the
  * queue. That is, in reverse order of which they were added, in order
  * to work properly with ResampleFilter.
  *
  * Imagine a user wants to apply an Affine transform followed by a Deformation
  * Field (DF) transform. He adds the Affine, then the DF. Because the user
  * typically conceptualizes a transformation as being applied from the Moving
  * image to the Fixed image, this makes intuitive sense. But since the
  * ResampleFilter expects to transform from the Fixed image to the Moving
  * image, the transforms are applied in reverse order of addition, i.e. from
  * the back of the queue, and thus, DF then Affine.
  */
  virtual OutputPointType TransformPoint( const InputPointType & inputPoint ) const{OutputPointType o; return o;}; 

  /**  Method to transform a vector. */
  // using Superclass::TransformVector;
  virtual OutputVectorType TransformVector(const InputVectorType &) const{OutputVectorType o; return o;};

  // virtual OutputVnlVectorType TransformVector(const InputVnlVectorType & inputVector) const;

  // virtual OutputVectorPixelType TransformVector(const InputVectorPixelType & inputVector ) const;

  // virtual OutputVectorType TransformVector(const InputVectorType & inputVector,
  //                                          const InputPointType & inputPoint ) const;

  // virtual OutputVnlVectorType TransformVector(const InputVnlVectorType & inputVector,
  //                                             const InputPointType & inputPoint ) const;

  // virtual OutputVectorPixelType TransformVector(const InputVectorPixelType & inputVector,
  //                                               const InputPointType & inputPoint ) const;

  // /**  Method to transform a CovariantVector. */
  // using Superclass::TransformCovariantVector;
  // virtual OutputCovariantVectorType TransformCovariantVector(const InputCovariantVectorType &) const;

  // virtual OutputVectorPixelType TransformCovariantVector(const InputVectorPixelType &) const;

  // virtual OutputCovariantVectorType TransformCovariantVector(const InputCovariantVectorType & inputVector,
  //                                                            const InputPointType & inputPoint ) const;

  // virtual OutputVectorPixelType TransformCovariantVector(const InputVectorPixelType & inputVector,
  //                                                        const InputPointType & inputPoint ) const;

  // /** Method to transform a DiffusionTensor3D */
  // using Superclass::TransformDiffusionTensor3D;
  // virtual OutputDiffusionTensor3DType TransformDiffusionTensor3D(
  //   const InputDiffusionTensor3DType & inputTensor) const;

  // virtual OutputVectorPixelType TransformDiffusionTensor3D(
  //   const InputVectorPixelType & inputTensor) const;

  // virtual OutputDiffusionTensor3DType TransformDiffusionTensor3D(
  //   const InputDiffusionTensor3DType & inputTensor,
  //   const InputPointType & inputPoint ) const;

  // virtual OutputVectorPixelType TransformDiffusionTensor3D(
  //   const InputVectorPixelType & inputTensor,
  //   const InputPointType & inputPoint ) const;

  // /** Method to transform a SymmetricSecondRankTensor */
  // using Superclass::TransformSymmetricSecondRankTensor;
  // virtual OutputSymmetricSecondRankTensorType TransformSymmetricSecondRankTensor(
  //   const InputSymmetricSecondRankTensorType & inputTensor) const;

  // virtual OutputVectorPixelType TransformSymmetricSecondRankTensor(
  //   const InputVectorPixelType & inputTensor) const;

  // virtual OutputSymmetricSecondRankTensorType TransformSymmetricSecondRankTensor(
  //   const InputSymmetricSecondRankTensorType & inputTensor,
  //   const InputPointType & inputPoint ) const;

  // virtual OutputVectorPixelType TransformSymmetricSecondRankTensor(
  //   const InputVectorPixelType & inputTensor,
  //   const InputPointType & inputPoint ) const;

  /**
   * Compute the Jacobian with respect to the parameters for the compositie
   * transform using Jacobian rule. See comments in the implementation.
   */
  virtual void ComputeJacobianWithRespectToParameters(const InputPointType  & p, JacobianType & j) const {};
  virtual void ComputeJacobianWithRespectToPosition(const InputPointType &, JacobianType &) const {};



  /** Specific Transform Types for PolyAffine **/
  typedef ConstantVelocityFieldTransform<TScalar, NDimensions>
    ConstantVelocityFieldTransformType;
  typedef MatrixOffsetTransformBase<TScalar, NDimensions, NDimensions> MatrixOffsetTransformBaseType;
  typedef typename ConstantVelocityFieldTransformType::ConstantVelocityFieldType ConstantVelocityFieldType;
  typedef typename ConstantVelocityFieldType::PixelType VelocityPixelType;
  typedef typename ConstantVelocityFieldTransformType::DisplacementFieldType DisplacementFieldType;



  /** Trajectory Collision **/
  // unsigned int m_SampleDistance;
  typedef std::vector<InputPointType> PointVectorType;
  typedef std::vector<PointVectorType> VectorOfPointVectorType;
  typedef std::vector<unsigned int> PointLabelVectorType; 
  // typedef std::vector<VectorOfPointVectorType> VectorOfVectorOfPointVectorType;
  double m_InterTrajectoryDistThres;
  itkSetMacro(InterTrajectoryDistThres, double);
  itkGetMacro(InterTrajectoryDistThres, double);

  PointVectorType m_ControlPointList;
  PointLabelVectorType m_ControlPointLabelList;

  unsigned int m_NumberOfTimeSteps;
  itkSetMacro(NumberOfTimeSteps, unsigned int);
  itkGetMacro(NumberOfTimeSteps, unsigned int);

  std::vector<float> m_CollisionTimeList;
  // points[pts_idx, time]
  VectorOfPointVectorType m_ControlPointListTrajectory;
  
  double m_BoundarySigma;
  itkSetMacro(BoundarySigma, double);
  itkGetMacro(BoundarySigma, double);




  /** Individual Affine Transform **/
  // std::vector<typename MatrixOffsetTransformBaseType::Pointer>  
  //   LocalAffineList;


  typedef GSLocalAffineTransform<TScalar, NDimensions> LocalAffineTransformType;

  std::vector<typename LocalAffineTransformType::Pointer> m_LocalAffineList;



  typedef unsigned int LabelType;
  typedef itk::Image<LabelType, NDimensions> LabelImageType;
  typedef typename LabelImageType::Pointer LabelImagePointerType;
  typename LabelImageType::Pointer m_LabelImage;

  /** Final Stationary Velocity Transform **/
  // typedef typename ConstantVelocityFieldTransformType::Pointer m_FieldTransfrom;

  typedef typename LabelImageType::IndexType SampleSpacingType;

  SampleSpacingType m_SampleSpacing;


  typedef typename ConstantVelocityFieldTransformType::Pointer ConstantVelocityFieldTransformPointerType;


  typedef Image<unsigned char, NDimensions> BinaryImageType;
  typedef typename BinaryImageType::Pointer BinaryImagePointerType;
  typedef Image<float, NDimensions> FloatImageType;
  typedef typename FloatImageType::Pointer FloatImagePointerType;
  typedef Image<float, NDimensions> DistanceImageType;
  typedef typename DistanceImageType::Pointer DistanceImagePointerType;

  struct BandPointSample{
    typename LabelImageType::IndexType index; // point index in the band
    LabelType label; // label of the point in the voronoi, s.t neighbor_weights should be dominated by the affine of this label
    std::vector<LabelType> neighbor_labels; // labels for nearest neighbors of the point of each affine mask inside the band radius
    std::vector<float> neighbor_distances;
    std::vector<float> neighbor_weights; // computed from neighbor_distances
  };

  typedef std::vector<BandPointSample> BandPointSampleVector;

  struct perCollision {
    LabelImagePointerType trajectoryImage;
    DistanceImagePointerType trajectoryDistanceImage;
    LabelImagePointerType trajectoryDistanceVoronoiImage;
    BinaryImagePointerType voronoiBandImage;
    BinaryImagePointerType voronoiEdgeImage;
    BandPointSampleVector bandPointSampleVector;
    ConstantVelocityFieldTransformPointerType velocityFieldTransform;
  };

  std::vector<perCollision> m_C;

//  std::vector<LabelImagePointerType> m_TrajectoryImageList;
//  std::vector<DistanceImagePointerType> m_TrajectoryDistanceImageList;
//  std::vector<LabelImagePointerType> m_TrajectoryDistanceVoronoiImageList;
//
//  std::vector<ConstantVelocityFieldTransformPointerType> m_VelocityFieldList;


  typedef ShapedNeighborhoodIterator<LabelImageType> SNIteratorType;


  typedef CompositeTransform<TScalar, NDimensions> CompositeTransformType;
  typedef typename CompositeTransformType::Pointer CompositeTransformPointerType;

  void ComputeFieldTransfromFromLocalAffineTransform();
  void SampleFromLabelMask(const LabelImagePointerType &, const SampleSpacingType &, PointVectorType &, PointLabelVectorType &);
  void ConvertAffineTransformsToLogarithm();
  void PrecomputeTrajectoryAndCheckCollision();
  void PrecomputeTrajectoryAndCheckCollision2(); // use rendering each point to a ball to check collision
  bool AddSphereToLabelMask(LabelImagePointerType &collision_image, const InputPointType &pt, const unsigned int radius, const double inter_trajectory_dist2_thres, const unsigned int labelK);
  void GetStationaryVieldCopyPasteDecreasing(unsigned int collision_time, ConstantVelocityFieldTransformPointerType &vfield);
  void ScatterLabeledPointsIntoTrajectoryImage(unsigned int collision_time, LabelImagePointerType &traj, SNIteratorType &it);

  // bool AddSphereToLabelMaskIgnoreCollision(LabelImagePointerType &collision_image, const IndexType &pt, const unsigned int radius, const unsigned int labelK){
  // SNIteratorType CreateBallShapedNeighborhoodIterator(LabelImagePointerType &collision_image, const unsigned int radius);
  void FillBallShapedNeighborhoodIterator(const unsigned int radius, SNIteratorType &it);
  void GetTrajectoryVoronoiBand(const LabelImagePointerType &voronoi, unsigned int voronoi_band_radius, BinaryImagePointerType &dL, BinaryImagePointerType &edge);
  void AssignWeightsToPointsInsideBand(const BinaryImagePointerType &band, const BinaryImagePointerType &edge, unsigned int collision_time);
  void AssignWeightsToPointsInsideBand2(const BinaryImagePointerType &band, const BinaryImagePointerType &edge, unsigned int collision_time);
  void MixingVelocityFieldsFromEveryLabel(unsigned int collision_time);
  void SuppressVelocityCloseToBoundary(unsigned int collision_time);

  // m_CompositeTransform is the concatenation of all the velocity field transform
  // note that: pay attention to the order when adding to composite transform
  // the reverse order in the collision time?
  CompositeTransformPointerType m_CompositeTransform;
  void ComputeCompositeTransform();

  void ComputeTrajectoryOfAllTimeInOneImage();

  LabelImagePointerType m_TrajectoryOfAllTime;

  void ComputeWeightByDistanceToTrajectoryOfAllTime();
  DistanceImagePointerType m_TrajectoryOfAllTimeDistanceMap;
  FloatImagePointerType m_TrajectoryOfAllTimeWeightImage;

  double m_TrajectoryOfAllTimeDistanceSigma;
  itkSetMacro(TrajectoryOfAllTimeDistanceSigma, double);
  itkGetMacro(TrajectoryOfAllTimeDistanceSigma, double);

  unsigned int m_BandRadius; // voronoi_band_radius
  itkSetMacro(BandRadius, unsigned int);
  itkGetMacro(BandRadius, unsigned int);

  void SuppressVelocityCloseToTrajectory(unsigned int collision_time);

  void ScaleVelocityForExponential(unsigned int collision_time);

protected:

  GSPolyAffineTransform();
  virtual ~GSPolyAffineTransform(){};
  void PrintSelf( std::ostream& os, Indent indent ) const{};

  /** Clone the current transform */
  virtual typename LightObject::Pointer InternalClone() const{};


private:
  GSPolyAffineTransform( const Self & ); // purposely not implemented
  void operator=( const Self & );             // purposely not implemented

public:
  typedef typename LabelImageType::IndexType IndexType;
  float euclidean_distance_square(const IndexType &index1, const IndexType &index2)
  {
    float d2 = 0;
    for(unsigned int d=0; d<InputDimension; d++) d2 += (index1[d]-index2[d])*(index1[d]-index2[d]);
    return d2;
  }

};

// some my temp helper functions for generic image handling

template<class InputImagePointerType, class OutputImagePointerType>
void AllocateNewImageOfSameDimension(const InputImagePointerType &input, OutputImagePointerType &output)
{
  typedef typename OutputImagePointerType::ObjectType OutputImageType;

  std::cout << "input region:" << input->GetBufferedRegion() << std::endl;

  output = OutputImageType::New();
  output->CopyInformation( input );
  output->SetRequestedRegion( input->GetRequestedRegion() );
  output->SetBufferedRegion( input->GetBufferedRegion() );
  output->Allocate();
  output->FillBuffer(0);

  std::cout << "input region:" << input->GetBufferedRegion() << std::endl;
}

// use this when OutputImage is an vector Image of Vector
template<class InputImagePointerType, class OutputImagePointerType>
void AllocateNewVectorImageOfSameDimension(const InputImagePointerType &input, OutputImagePointerType &output)
{
  typedef typename OutputImagePointerType::ObjectType OutputImageType;

  std::cout << "input region:" << input->GetBufferedRegion() << std::endl;

  typedef typename OutputImageType::PixelType PixelType;
  PixelType zeroVector( itk::NumericTraits<PixelType>::ZeroValue() );
  output = OutputImageType::New();
  output->CopyInformation( input );
  output->SetRequestedRegion( input->GetRequestedRegion() );
  output->SetBufferedRegion( input->GetBufferedRegion() );
  output->Allocate();
  output->FillBuffer(zeroVector);

  std::cout << "input region:" << input->GetBufferedRegion() << std::endl;
}


template<class LabelImagePointerType, class DistanceImagePointerType>
void GetDistanceTransformFromLabelImage(const LabelImagePointerType &label, DistanceImagePointerType &distance, LabelImagePointerType &voronoi)
{

  typedef typename LabelImagePointerType::ObjectType LabelImageType;
  typedef typename DistanceImagePointerType::ObjectType DistanceImageType;
  typedef DanielssonDistanceMapImageFilter<LabelImageType, DistanceImageType > FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(label);
  filter->InputIsBinaryOff();

  std::cout << "in GetDistanceTransformFromLabelImage, Update()" << std::endl;
  filter->Update();

  std::cout << "get distance map out" << std::endl;
  distance = filter->GetOutput();
  std::cout << "get voronoi label out" << std::endl;
  voronoi = filter->GetVoronoiMap();
  std::cout << "out 3" << std::endl;


}

float compute_boundary_suppresion_weight(float d, float s)
{
  float w = 0;

  w = 1 - vcl_exp(-1 * (d/s));

  return w;
}

float compute_weight_from_distance_to_boundary(float pd, double th)
{
  double w = 0;

  if (pd <= -1 * th )
    w = 1.0;
  else if (pd < th)
    w = pd / (-2*th) + 0.5;
  else
    w = 0;

  return w;
}

template<class ImageType, class String>
void GSWriteImageRawPointer(const ImageType *img, const String &str)
{
  typedef ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(img);
  writer->SetFileName(str);
  writer->Update();
}

template<class ImagePointerType, class String>
void GSWriteImage(const ImagePointerType &img, const String &str)
{
  typedef typename ImagePointerType::ObjectType ImageType;
  typedef ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(img);
  writer->SetFileName(str);
  writer->Update();
}




} // end namespace itk

#if ITK_TEMPLATE_TXX
#include "itkGSPolyAffineTransform.hxx"
#endif

#endif // __itkGSPolyAffineTransform_h
