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

template <unsigned int Dimension>
int itkPolyAffineTransformTest(int argc, char *argv[])
{
  char *fileName = argv[2];

  PicslImageHelper::m_FilePath[0] = '\0';
  if (argc >= 4)
    {
    strcpy(PicslImageHelper::m_FilePath, argv[3]);
    }

  PicslImageHelper::m_Debug = false;
  if (argc >= 5 && strcmp(argv[4], "true")==0)
    {
    PicslImageHelper::m_Debug = true;
    }

  //create a deformation field transform
  //typedef TranslationTransform<double, Dimension>

  typedef typename itk::PolyAffineTransform<double, Dimension> PolyAffineTransformType;
  typedef typename PolyAffineTransformType::LocalAffineTransformType LocalAffineTransformType;
  typedef typename PolyAffineTransformType::DisplacementFieldType DisplacementFieldType;
  typedef typename LocalAffineTransformType::MaskImageType MaskImageType;

  typedef typename itk::GroupSpatialObject< Dimension > SceneType;
  typedef typename itk::BoxSpatialObject< Dimension > BoxType;

  typename SceneType::Pointer scene1 = SceneType::New();
  typename SceneType::Pointer scene2 = SceneType::New();
  typename BoxType::Pointer box1 = BoxType::New();
  typename BoxType::Pointer box2 = BoxType::New();
  scene1->AddSpatialObject(box1);
  scene2->AddSpatialObject(box2);
  
  typename BoxType::SizeType boxsize1;
  typename BoxType::SizeType boxsize2;
  boxsize1.Fill(20);
  box1->SetSize( boxsize1 );
  boxsize2.Fill(20);
  box2->SetSize( boxsize2 );

  typename BoxType::TransformType::OffsetType boxOffset1;
  typename BoxType::TransformType::OffsetType boxOffset2;
  //not overlapping below
  //double offset1[] = {50.0, 10.0, 50.0};
  //double offset2[] = {10.0, 50.0, 50.0};
  //crossing below
  //double offset1[] = {60.0, 20.0, 50.0};
  //double offset2[] = {20.0, 60.0, 50.0};
  //overlapping below
  //double offset1[] = {10.0, 100.0, 50.0};
  //double offset2[] = {100.0, 100.0, 50.0};
  //from alex below
  double offset1[] = {80-60.0, 30-10.0, 20.0};
  double offset2[] = {80+40.0, 30-10.0, 120.0};
  for (unsigned int d=0; d<Dimension; d++)
    {
    boxOffset1[d] = offset1[d];
    boxOffset2[d] = offset2[d];
    }
  box1->GetObjectToParentTransform()->SetOffset( boxOffset1 );
  box1->ComputeObjectToWorldTransform();
  box2->GetObjectToParentTransform()->SetOffset( boxOffset2 );
  box2->ComputeObjectToWorldTransform();

  typename MaskImageType::SizeType size;
  size.Fill(160);

  typename PolyAffineTransformType::Pointer polyTransform = PolyAffineTransformType::New();
  typename LocalAffineTransformType::Pointer localTransform1 = LocalAffineTransformType::New();
  typename LocalAffineTransformType::Pointer localTransform2 = LocalAffineTransformType::New();

  typedef itk::Matrix<double, Dimension, Dimension> MatrixType;
  typedef itk::Vector<double, Dimension> VectorType;
  VectorType affineOffset1, affineOffset2;
  //not overlapping below
  //double translation1[] = {50, 0, 30};
  //double translation2[] = {0, 50, 30};
  //crossing below
  //double translation1[] = {0, 80, 30};
  //double translation2[] = {80, 0, 30};
  //overlapping below
  //double translation1[] = {100, 0, 30};
  //double translation2[] = {0, -80, 30};
  //from alex below
  double translation1[] = {100, 0, 20};
  double translation2[] = {0, 50, -20};
  for (unsigned int d=0; d<Dimension; d++)
    { 
    affineOffset1[d] = translation1[d];
    affineOffset2[d] = translation2[d];
    }
  localTransform1->SetOffset(affineOffset1);
  localTransform2->SetOffset(affineOffset2);

  localTransform1->template ComputeFixedMaskImageFromSpatialObject< SceneType >(scene1, size);
  localTransform2->template ComputeFixedMaskImageFromSpatialObject< SceneType >(scene2, size);

  polyTransform->AddLocalAffineTransform(localTransform1);
  polyTransform->AddLocalAffineTransform(localTransform2);
  polyTransform->SetTimeStampLog(8);

  typename DisplacementFieldType::Pointer displacementField = polyTransform->GetDisplacementField();

  PicslImageHelper::WriteResultImage<DisplacementFieldType>(displacementField, fileName);
  PicslImageHelper::WriteResultDisplacementField<DisplacementFieldType>(displacementField, fileName);

  std::cout << "Test PASSED." << std::endl;

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " Dimension=2|3 OutputDisplacementFieldFile [OutputDirectory] [DebugFiles=true|false]" << std::endl;
    std::cerr << "Example: " << argv[0];
    std::cerr << " 2 outfield.nii outdir true" << std::endl;
    return EXIT_FAILURE;
    }

  int Dimension = atoi(argv[1]);
  if (Dimension == 2)
    {
    itkPolyAffineTransformTest<2>(argc, argv);
    }
  else if (Dimension == 3)
    {
    itkPolyAffineTransformTest<3>(argc, argv);
    }
  else
    {
    std::cerr << "Dimension " << Dimension << " is not supported." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
