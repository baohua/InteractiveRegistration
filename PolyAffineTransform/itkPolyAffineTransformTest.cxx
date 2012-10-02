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
/**
 * Test program for DemonImageToImageObjectMetric and
 * GradientDescentObjectOptimizer classes.
 *
 * Tests are disabled since it requires some image files as input.
 */

#include <string.h>

#include "itkPolyAffineTransform.h"
#include "itkImageFileWriter.h"
#include "itkGroupSpatialObject.h"
#include "itkBoxSpatialObject.h"
#include "itkPicslImageHelper.h"

using namespace itk;

int itkPolyAffineTransformTest(int argc, char *argv[])
{
  if( argc < 2 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " OutputDisplacementFieldFile ";
    return EXIT_FAILURE;
    }

  //create a deformation field transform
  //typedef TranslationTransform<double, Dimension>
  const int Dimension = 2;

  typedef itk::PolyAffineTransform<double, Dimension> PolyAffineTransformType;
  typedef PolyAffineTransformType::LocalAffineTransformType LocalAffineTransformType;
  typedef PolyAffineTransformType::DisplacementFieldType DisplacementFieldType;
  typedef LocalAffineTransformType::MaskImageType MaskImageType;
  typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;

  PolyAffineTransformType::Pointer polyTransform = PolyAffineTransformType::New();
  LocalAffineTransformType::Pointer localTransform1 = LocalAffineTransformType::New();
  LocalAffineTransformType::Pointer localTransform2 = LocalAffineTransformType::New();

  typedef  itk::Matrix<double, Dimension, Dimension> MatrixType;
  typedef  itk::Vector<double, Dimension> VectorType;
  VectorType affineOffset1, affineOffset2;
  double translation1[] = {30, 30, 30};
  double translation2[] = { 0, 30,  0};
  for (unsigned int d=0; d<Dimension; d++)
    {
    affineOffset1[d] = translation1[d];
    affineOffset2[d] = translation2[d];
    }
  localTransform1->SetOffset(affineOffset1);
  localTransform2->SetOffset(affineOffset2);

  typedef itk::GroupSpatialObject< Dimension > SceneType;
  typedef itk::BoxSpatialObject< Dimension > BoxType;

  SceneType::Pointer scene1 = SceneType::New();
  SceneType::Pointer scene2 = SceneType::New();
  BoxType::Pointer box1 = BoxType::New();
  BoxType::Pointer box2 = BoxType::New();
  scene1->AddSpatialObject(box1);
  scene2->AddSpatialObject(box2);
  
  BoxType::SizeType boxsize1;
  BoxType::SizeType boxsize2;
  boxsize1.Fill(30);
  box1->SetSize( boxsize1 );
  boxsize2.Fill(30);
  box2->SetSize( boxsize2 );

  BoxType::TransformType::OffsetType boxOffset1;
  BoxType::TransformType::OffsetType boxOffset2;
  double offset1[] = {50.0, 10.0, 10.0};
  double offset2[] = {10.0, 50.0, 50.0};
  for (unsigned int d=0; d<Dimension; d++)
    {
    boxOffset1[d] = offset1[d];
    boxOffset2[d] = offset2[d];
    }
  box1->GetObjectToParentTransform()->SetOffset( boxOffset1 );
  box1->ComputeObjectToWorldTransform();
  box2->GetObjectToParentTransform()->SetOffset( boxOffset2 );
  box2->ComputeObjectToWorldTransform();

  MaskImageType::SizeType size;
  size.Fill(128);
  localTransform1->ComputeFixedMaskImageFromSpatialObject<SceneType>(scene1, size);
  localTransform1->ComputeMovingMaskImageAndDenseFixedPointSet();
  itk::PicslImageHelper::WriteImage<LocalAffineTransformType::MaskImageType>(localTransform1->GetFixedMaskImage(), "tmpFixedMask0.nii");
  itk::PicslImageHelper::WriteImage<LocalAffineTransformType::MaskImageType>(localTransform1->GetMovingMaskImage(), "tmpMovingMask0.nii");

  localTransform2->ComputeFixedMaskImageFromSpatialObject<SceneType>(scene2, size);
  localTransform2->ComputeMovingMaskImageAndDenseFixedPointSet();
  itk::PicslImageHelper::WriteImage<LocalAffineTransformType::MaskImageType>(localTransform2->GetFixedMaskImage(), "tmpFixedMask1.nii");
  itk::PicslImageHelper::WriteImage<LocalAffineTransformType::MaskImageType>(localTransform2->GetMovingMaskImage(), "tmpMovingMask1.nii");

  polyTransform->AddLocalAffineTransform(localTransform1);
  polyTransform->AddLocalAffineTransform(localTransform2);
  polyTransform->SetTimeStampLog(8);

  DisplacementFieldType::Pointer displacementField = polyTransform->GetDisplacementField();

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(displacementField);
  writer->SetFileName( argv[ 1 ] );
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject & e )
    {
    std::cerr << "Exception detected while generating displacement field" << argv[1];
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
    }
  
  itk::PicslImageHelper::WriteDisplacementField<DisplacementFieldType>
    (displacementField, argv[1]);

  std::cout << "Test PASSED." << std::endl;

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  itkPolyAffineTransformTest(argc, argv);
}
