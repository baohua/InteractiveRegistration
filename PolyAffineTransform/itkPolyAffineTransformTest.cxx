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

#include "itkPolyAffineTransform.h"
#include "itkImageFileWriter.h"
#include "itkGroupSpatialObject.h"
#include "itkBoxSpatialObject.h"

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
  affineOffset1[0] = -10;
  affineOffset1[1] = 0;
  affineOffset2[0] = 10;
  affineOffset2[1] = 0;
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
  
  BoxType::SizeType  boxsize1;
  BoxType::SizeType  boxsize2;
  boxsize1[0] = 30;
  boxsize1[1] = 30;
  box1->SetSize( boxsize1 );
  boxsize2[0] = 30;
  boxsize2[1] = 30;
  box2->SetSize( boxsize2 );

  BoxType::TransformType::OffsetType offset1;
  BoxType::TransformType::OffsetType offset2;
  offset1[0] = 50.0;
  offset1[1] = 10.0;
  box1->GetObjectToParentTransform()->SetOffset( offset1 );
  box1->ComputeObjectToWorldTransform();
  offset2[0] = 10.0;
  offset2[1] = 50.0;
  box2->GetObjectToParentTransform()->SetOffset( offset2 );
  box2->ComputeObjectToWorldTransform();

  MaskImageType::SizeType size;
  size.Fill(100);
  localTransform1->SetMaskImageFromSpatialObject<SceneType>(scene1, size);
  localTransform2->SetMaskImageFromSpatialObject<SceneType>(scene2, size);

  polyTransform->AddLocalAffineTransform(localTransform1);
  polyTransform->AddLocalAffineTransform(localTransform2);
  polyTransform->SetTimeStampLog(8);

  DisplacementFieldType::Pointer displacementField = polyTransform->GetExponentialDisplacementField();

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

  std::cout << "Test PASSED." << std::endl;

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  itkPolyAffineTransformTest(argc, argv);
}
