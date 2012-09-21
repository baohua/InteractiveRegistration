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
{}

/** Constructor with default arguments */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::LocalAffineTransform(unsigned int parametersDimension):
  Superclass(parametersDimension)
{}

/** Constructor with explicit arguments */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::LocalAffineTransform(const MatrixType & matrix,
                                                             const OutputVectorType & offset):
  Superclass(matrix, offset)
{}

/**  Destructor */
template< class TScalarType, unsigned int NDimensions >
LocalAffineTransform< TScalarType, NDimensions >::
~LocalAffineTransform()
{
  return;
}

/** Compute principal logorithm of the homegeneous matrix */
template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >::ComputePrincipalLogorithm()
{
  this->m_VelocityAffineTransform = VelocityAffineTransformType::New();

  vnl_matrix<TScalarType> homoMat(NDimensions+1, NDimensions+1);
  MatrixType mat = this->GetMatrix();
  OutputVectorType offset = GetOffset();

  for (unsigned int i=0; i<NDimensions; i++)
    {
    for (unsigned int j=0; j<NDimensions; j++)
      {
      homoMat[i][j] = mat[i][j];
      }
    homoMat[i][NDimensions] = offset[i];
    }
  homoMat[NDimensions][NDimensions] = 1;

  vnl_matrix<TScalarType> logMat = sdtools::GetLogarithm(homoMat);

  MatrixType vmat;
  OutputVectorType voffset;
  for (unsigned int i=0; i<NDimensions; i++)
    {
    for (unsigned int j=0; j<NDimensions; j++)
      {
        vmat[i][j] = logMat[i][j];
      }
    voffset[i] = logMat[i][NDimensions];
    }

  this->m_VelocityAffineTransform->SetMatrix(vmat);
  this->m_VelocityAffineTransform->SetOffset(voffset);

}

template< class TScalarType, unsigned int NDimensions >
typename LocalAffineTransform< TScalarType, NDimensions >::VelocityAffineTransformType*
LocalAffineTransform< TScalarType, NDimensions >::GetVelocityAffineTransform()
{
  if (this->m_VelocityAffineTransform.IsNull())
    {
    this->ComputePrincipalLogorithm();
    }
  return this->m_VelocityAffineTransform.GetPointer();
}

template< class TScalarType, unsigned int NDimensions >
void
LocalAffineTransform< TScalarType, NDimensions >::AddControlPoint(InputPointType point)
{
  if (this->m_ControlPointSet.IsNull())
    {
    this->m_ControlPointSet = PointSet::New();
    }
  this->m_ControlPointSet.push_back(point);
}

template< class TScalarType, unsigned int NDimensions >
template< class TSpatialObject > 
void 
LocalAffineTransform< TScalarType, NDimensions >
::SetMaskImageFromSpatialObject(TSpatialObject *spatialObject,
                                SizeType size)
{
  typedef itk::SpatialObjectToImageFilter<TSpatialObject,MaskImageType> SpatialObjectToImageFilterType;
  SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
  imageFilter->SetInput(spatialObject);
  imageFilter->SetInsideValue(1);
  imageFilter->SetOutsideValue(0);

  imageFilter->SetSize(size);

  imageFilter->Update();
  this->m_MaskImage = imageFilter->GetOutput();   
}
template< class TScalarType, unsigned int NDimensions >
template< class TSpatialObject > 
void 
LocalAffineTransform< TScalarType, NDimensions >
::SetMaskImageFromSpatialObject(TSpatialObject *spatialObject,
                                SizeType size,
                                PointType origin,
                                DirectionType direction,
                                SpacingType spacing)
{
  typedef itk::SpatialObjectToImageFilter<TSpatialObject,MaskImageType> SpatialObjectToImageFilterType;
  SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
  imageFilter->SetInput(spatialObject);
  imageFilter->SetInsideValue(1);
  imageFilter->SetOutsideValue(0);

  imageFilter->SetOrigin(origin);
  imageFilter->SetSpacing(spacing);
  imageFilter->SetDirection(direction);
  imageFilter->SetSize(size);

  imageFilter->Update();
  this->m_MaskImage = imageFilter->GetOutput();
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
