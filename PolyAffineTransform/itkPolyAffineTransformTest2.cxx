/*
 * itkPolyAffineTransformTest2.cxx
 *
 *  Created on: Nov 12, 2012
 *      Author: songgang
 */

/**
 * My test script for Poly-affine transform
 *
 * 1. Input test: two ways to define regions
 *  1) Use a image label mask ==> convert to a point set with a sampling radius
 *  2) Use a point set with labels
 *  3) Use spatial object??
 *
 *  Store: an image label (for distance transform etc) + point set (for
 *  iteration)
 *
 * 2. Input test: use variant sub-types of affine transform as local affine
 * transform
 *
 * 3. Output test:
 *  1) Generate the weight map from the local region mask
 *  1) Generate the velocity field
 *  2) Generate the displacement
 *
 * 4. Gradient test: need to multi-thread, take advantage of affine
 * transform class
 *  Add bridge from local metric change to affine trnasform (similar to
 *  composite transform chain rule)
 *  1) Compare the gradient test of single affine transform with using the
 *  full image as the mask
 *
 * 5. Apple the transform to an image;
 *
 */


#include "itkCompositeTransform.h"
#include "itkImageMaskSpatialObject.h"

#include "itkAffineTransform.h"

using namespace itk;

template<int Dimension, class TMaskImagePointerType>
void GetTestMaskImage(TMaskImagePointerType &mask)
{
  typedef itk::ImageMaskSpatialObject<Dimension> ImageMaskSpatialObjectType;
  typedef typename ImageMaskSpatialObjectType::PixelType MaskPixelType;
  typedef typename ImageMaskSpatialObjectType::ImageType MaskImageType;

//  typename MaskImageType::Pointer mask = MaskImageType::New();
  typename MaskImageType::SizeType size;
  size.Fill(50);
  typename MaskImageType::IndexType index;
  index.Fill(0);
  typename MaskImageType::RegionType region;

  region.SetSize(size);
  region.SetIndex(index);

  mask->SetRegions(region);
  mask->Allocate();

  MaskPixelType p = itk::NumericTraits< MaskPixelType >::Zero;

  mask->FillBuffer(p);

  MaskImageType::RegionType insideRegion;
  MaskImageType::SizeType insideSize = {{ 30, 30, 30 }};
  MaskImageType::IndexType insideIndex = {{ 10, 10, 10 }};
  insideRegion.SetSize( insideSize );
  insideRegion.SetIndex( insideIndex );

  Iterator it( image, insideRegion );
  it.GoToBegin();

  while ( !it.IsAtEnd() ) {
    it.Set( 1 );
    ++it;
  }

}

void ComputeDeformationField()
{
  //try to duplicate the matlab code


}

template <unsigned int Dimension>
int itkPolyAffineTransformTest2(int argc, char *argv[])
{
  typedef itk::ImageMaskSpatialObject<Dimension> ImageMaskSpatialObjectType;
  typedef typename ImageMaskSpatialObjectType::ImageType MaskImageType;

  typename MaskImageType::Pointer mask = MaskImageType::New();

  GetTestMaskImage(mask);

  typename ImageMaskSpatialObjectType::Pointer maskSO = ImageMaskSpatialObjectType::New();
  maskSO->SetImage(mask);

  // load the affine transform
  typedef itk::AffineTransform  AffineTransformType;
  AffineTransformType::Pointer aff1 = AffineTransformType::New();
  AffineTransformType::Pointer aff2 = AffineTransformType::New();

  GetTestLeftAffineTransform(aff1);
  GetTestRightAffineTransform(aff2);

  typedef itk::MatPolyAffineTransform PolyAffineTransformType;
  PolyAffineTransformType::Pointer paff = PolyAffineTransformType::New();

  paff->SetRegionLabel(mask);
  unsigned int label = 1;
  paff->SetLocalTransform(label, aff1);
  label = 2;
  paff->SetLocalTransform(label, aff2);

  paff->ComputeDeformationField();

  // for debug use
  paff->GetWeightImage(1);
  paff->GetWeightImage(2);
  paff->GetVelocityField();

  paff->GetJacobianWithRespectToParameter();
  paff->GetJacobianWithRespectToPosition()

  // add the registration pipeline;
  metric->UsePointSetSampling();
  registration->SetMetric(metric);
  registration->SetTransform(paff);











}

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
    std::cerr << " Dimension 2|3 InputLabelMask AffineFileForLabel_1 AffineFileFor Lable_2 ... (to K, K=maximum of InputLabelMask)"
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



