#include "itkGSPolyAffineTransform.h"
#include "itkImageFileReader.h"

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
    offset[0] = 12; offset[1] = 20; offset[2] = 30;
  }
  else {

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
    matrix(0,0) = 1; matrix(0,1) = 1; matrix(0,2) = 0;
    matrix(1,0) = 0; matrix(1,1) = 1; matrix(1,2) = 0;
    matrix(2,0) = 0; matrix(2,1) = 0; matrix(2,2) = 1;
    offset[0] = -12; offset[1] = -20; offset[2] = 25;
  }
  else {

  }

  aff->SetMatrix(matrix);
  aff->SetOffset(offset);
}


template <unsigned int Dimension>
int itkGSPolyAffineTransformTest(int argc, char **argv)
{
	typedef typename itk::GSPolyAffineTransform<double, Dimension>
	PolyAffineTransformType;

	typename PolyAffineTransformType::Pointer paff = PolyAffineTransformType::New();


	typedef typename PolyAffineTransformType::LabelImageType LabelImageType;

	typename itk::ImageFileReader<LabelImageType>::Pointer reader = itk::ImageFileReader<LabelImageType>::New();
	reader->SetFileName(argv[1]);
	reader->Update();

	paff->m_LabelImage = reader->GetOutput();

	typename PolyAffineTransformType::SampleSpacingType sample_spacing;
	sample_spacing[0] = 2;
	sample_spacing[1] = 2;
	sample_spacing[2] = 2;

	paff->m_SampleSpacing = sample_spacing;

	std::cout << "imgsz" << paff->m_LabelImage->GetLargestPossibleRegion().GetSize() << std::endl;

  // load the affine transform
  typedef itk::AffineTransform<double, 3>  AffineTransformType;
  AffineTransformType::Pointer aff1 = AffineTransformType::New();
  AffineTransformType::Pointer aff2 = AffineTransformType::New();

  GetTestLeftAffineTransform(aff1);
  GetTestRightAffineTransform(aff2);

  paff->AddTransform(aff1);
  paff->AddTransform(aff2);

  paff->SetInterTrajectoryDistThres(5.2);
  paff->SetNumberOfTimeSteps(50);

	paff->ComputeFieldTransfromFromLocalAffineTransform();

}



int main(int argc, char **argv)
{
	itkGSPolyAffineTransformTest<3>(argc, argv);
	return EXIT_SUCCESS;
}
