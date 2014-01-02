#include "itkGSPolyAffineTransform.h"
#include "itkImageFileReader.h"
#include "itkPicslImageHelper.h"

template<unsigned int T>
struct DimTypeHelper
{
};

template<class ImageTypePointer >
void GetTestLabelImage(ImageTypePointer &img, DimTypeHelper<2>)
{
  typedef typename ImageTypePointer::ObjectType ImageType;
  img = ImageType::New();
  typename ImageType::RegionType region;
  typename ImageType::IndexType index;
  typename ImageType::SizeType size;
  size[0] = 301;
  size[1] = 301;
  index[0] = -150;
  index[1] = -150;
  region.SetSize(size);
  region.SetIndex(index);

  img->SetRegions(region);
  img->Allocate();

  itk::ImageRegionIteratorWithIndex<ImageType> it(img, img->GetBufferedRegion());
  for(it.GoToBegin(); ! it.IsAtEnd(); ++it) {
    typename ImageType::IndexType index1 = it.GetIndex();
    typename ImageType::PixelType label = 0;

    // assign label to two rectangle region
    if (index1[0]>=-60 && index1[0]<=-40 && index1[1]>=-10 && index1[1]<=10 ) {
      label = 1;
    }
    else if (index1[0]>=40 && index1[0]<=60 && index1[1]>=-10 && index1[1]<=10 ) {
      label = 2;
    }

    // std::cout << "index=" << index1 << " label=" << label << std::endl;
    it.Set(label);
  }

}


template<class ImageTypePointer >
void GetTestLabelImage(ImageTypePointer &img, DimTypeHelper<3>)
{
  typedef typename ImageTypePointer::ObjectType ImageType;
  img = ImageType::New();
  typename ImageType::RegionType region;
  typename ImageType::IndexType index;
  typename ImageType::SizeType size;
  size[0] = 301;
  size[1] = 301;
  size[2] = 101;
  index[0] = -150;
  index[1] = -150;
  index[2] = -50;
  region.SetSize(size);
  region.SetIndex(index);

  img->SetRegions(region);
  img->Allocate();

  itk::ImageRegionIteratorWithIndex<ImageType> it(img, img->GetBufferedRegion());
  for(it.GoToBegin(); ! it.IsAtEnd(); ++it) {
    typename ImageType::IndexType index1 = it.GetIndex();
    typename ImageType::PixelType label = 0;

    // assign label to two rectangle region
    if (index1[0]>=-60 && index1[0]<=-40 && index1[1]>=-10 && index1[1]<=10 && index1[2]>=-10 && index1[2]<=10) {
      label = 1;
    }
    else if (index1[0]>=40 && index1[0]<=60 && index1[1]>=-10 && index1[1]<=10  && index1[2]>=-10 && index1[2]<=10) {
      label = 2;
    }

    // std::cout << "index=" << index1 << " label=" << label << std::endl;
    it.Set(label);
  }

}


template<class TAffineTransformPointer>
void GetTestLeftAffineTransform(TAffineTransformPointer &aff) {
  typedef typename TAffineTransformPointer::ObjectType AffineTransformType;
  const int Dimension = AffineTransformType::InputSpaceDimension;
  typename AffineTransformType::MatrixType matrix;
  typename AffineTransformType::OffsetType offset;

  if (Dimension == 3) {
    matrix(0,0) = 1; matrix(0,1) = 0; matrix(0,2) = 0;
    matrix(1,0) = 0; matrix(1,1) = 1; matrix(1,2) = 0;
    matrix(2,0) = 0; matrix(2,1) = 0; matrix(2,2) = 1;
    offset[0] = 100; offset[1] = 0; offset[2] = 0;
  }
  else { // Dimension == 2
    matrix(0,0) = 1; matrix(0,1) = 0;
    matrix(1,0) = 0; matrix(1,1) = 1;
    offset[0] = 100; offset[1] = 0;
    // offset[0] = 20; offset[1] = 0;
  }

  aff->SetMatrix(matrix);
  aff->SetOffset(offset);
}

template<class TAffineTransformPointer>
void GetTestRightAffineTransform(TAffineTransformPointer &aff) {
  typedef typename TAffineTransformPointer::ObjectType AffineTransformType;
  const int Dimension = AffineTransformType::InputSpaceDimension;
  typename AffineTransformType::MatrixType matrix;
  typename AffineTransformType::OffsetType offset;

  if (Dimension == 3) {
    matrix(0,0) = 1; matrix(0,1) = 0; matrix(0,2) = 0;
    matrix(1,0) = 0; matrix(1,1) = 1; matrix(1,2) = 0;
    matrix(2,0) = 0; matrix(2,1) = 0; matrix(2,2) = 1;
    offset[0] = 0; offset[1] = 50; offset[2] = 20;
  }
  else {// Dimension == 2
    matrix(0,0) = 1; matrix(0,1) = 0;
    matrix(1,0) = 0; matrix(1,1) = 1;
    offset[0] = 0; offset[1] = 50;
    // offset[0] = 0; offset[1] = 20;

  }

  aff->SetMatrix(matrix);
  aff->SetOffset(offset);
}


template <unsigned int Dimension>
int itkGSPolyAffineTransformTest(int argc, char **argv)
{
  typedef double TScalar;
	typedef typename itk::GSPolyAffineTransform<TScalar, Dimension>
	PolyAffineTransformType;

	typename PolyAffineTransformType::Pointer paff = PolyAffineTransformType::New();
	typedef typename PolyAffineTransformType::LabelImageType LabelImageType;
	typedef typename LabelImageType::Pointer LabelImagePointerType;

//	typename itk::ImageFileReader<LabelImageType>::Pointer reader = itk::ImageFileReader<LabelImageType>::New();
//	reader->SetFileName(argv[1]);
//	reader->Update();
//
//	paff->m_LabelImage = reader->GetOutput();

	const DimTypeHelper<Dimension> dimTypeHelper;
	GetTestLabelImage(paff->m_LabelImage, dimTypeHelper);
	itk::PicslImageHelper::WriteResultImage<LabelImageType>(paff->m_LabelImage, "label2d.mhd");


	typename PolyAffineTransformType::SampleSpacingType sample_spacing;
	for(unsigned int i=0; i<Dimension; i++) sample_spacing[i] = 2;

	paff->m_SampleSpacing = sample_spacing;

	std::cout << "imgsz" << paff->m_LabelImage->GetLargestPossibleRegion().GetSize() << std::endl;

  // load the affine transform
  typedef itk::AffineTransform<TScalar, Dimension>  AffineTransformType;
  typename AffineTransformType::Pointer aff1 = AffineTransformType::New();
  typename AffineTransformType::Pointer aff2 = AffineTransformType::New();

  GetTestLeftAffineTransform(aff1);
  GetTestRightAffineTransform(aff2);

  std::cout << "aff1: " << aff1 << std::endl;
  std::cout << "aff2: " << aff2 << std::endl;

  paff->AddTransform(aff1);
  paff->AddTransform(aff2);

  paff->SetInterTrajectoryDistThres(5.2);
  paff->SetNumberOfTimeSteps(50);
  paff->SetBoundarySigma(2);
  paff->SetTrajectoryOfAllTimeDistanceSigma(36);
  paff->SetBandRadius(6);

  std::cout << "start paff->ComputeFieldTransfromFromLocalAffineTransform();" << std::endl;

	paff->ComputeFieldTransfromFromLocalAffineTransform();
  return 0;
}

int main(int argc, char **argv)
{
	itkGSPolyAffineTransformTest<3>(argc, argv);
  // itkGSPolyAffineTransformTest<2>(argc, argv);
	return EXIT_SUCCESS;
}
