#ifndef __itkPicslImageHelper_hxx
#define __itkPicslImageHelper_hxx

#include "itkPicslImageHelper.h"
#include "itkImageFileWriter.h"

namespace itk
{

PicslImageHelper
::PicslImageHelper()
{
}

PicslImageHelper
::~PicslImageHelper()
{
}

void
PicslImageHelper::PrintSelf(std::ostream & os, Indent indent) const
{
}
template< class TRegion>
typename itk::VectorContainer<int, typename TRegion::IndexType>::Pointer 
PicslImageHelper
::GetCorners( const TRegion &region)
{
  const unsigned int NDimensions = TRegion::ImageDimension;
  int cornerNum = 1 << NDimensions;

  typedef itk::VectorContainer<int, typename TRegion::IndexType>
    VectorContainerType;
  typedef typename VectorContainerType::Pointer
    VectorContainerPointer;
  
  VectorContainerPointer cornerIndices = VectorContainerType::New();
  cornerIndices->Reserve(cornerNum);

  typename TRegion::SizeType size = region.GetSize();
  typename TRegion::IndexType firstCorner = region.GetIndex();

  for(int i=0; i<cornerNum; i++)
    {
    int bit;
    typename TRegion::IndexType corner;
    for (int d=0; d<NDimensions; d++)
      {
      bit = (int) (( i & (1 << d) ) != 0); // 0 or 1
      corner[d] = firstCorner[d] + bit * (size[d] - 1);
      }
    cornerIndices->SetElement(i, corner);
    }
  return cornerIndices;
}


template< class TImage>
void 
PicslImageHelper
::WriteImage(typename TImage::Pointer image, char *fname)
{
  if (!PicslImageHelper::m_Debug)
    {
    return;
    }
  typedef ImageFileWriter<TImage> WriterType;
  
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(image);
  writer->SetFileName( fname );

  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject & e )
    {
    std::cerr << "Exception detected while writing image " << fname;
    std::cerr << " : "  << e.GetDescription();
    }

  return;
}

template< class TField>
void 
PicslImageHelper
::WriteDisplacementField(typename TField::Pointer field, char *fname)
{
  if (!PicslImageHelper::m_Debug)
    {
    return;
    }
  const unsigned int ImageDimension = TField::ImageDimension;
  typedef typename TField::PixelType VectorType;
  typedef typename TField::IndexType IndexType;
  typedef typename VectorType::ValueType ScalarType;  

  typedef typename itk::Image<ScalarType, ImageDimension> ScalarImageType;
  typedef typename itk::ImageFileWriter<ScalarImageType> WriterType;

  char *dotPtr = strrchr(fname, '.');
  int  dotPos = dotPtr - fname;
  char prefix[256], splitName[256];

  strcpy(prefix, fname);
  prefix[dotPos] = '\0';

  typedef ImageRegionIteratorWithIndex< TField > IteratorType;
  IteratorType it( field, field->GetLargestPossibleRegion() );

  VectorType vec;
  IndexType index;
  
  std::vector<typename ScalarImageType::Pointer> imageVector(ImageDimension);
  for (unsigned int i=0; i<ImageDimension; i++)
    {
    typename ScalarImageType::Pointer image = ScalarImageType::New();
    imageVector[i] = image;
    imageVector[i]->CopyInformation(field);
    imageVector[i]->SetRegions(field->GetLargestPossibleRegion());
    imageVector[i]->Allocate();
    }

  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    vec = it.Get();
    index = it.GetIndex();
    for (unsigned int i=0; i<ImageDimension; i++)
      {
      imageVector[i]->SetPixel(index, vec[i]);
      }
    }

  for (unsigned int d=0; d<ImageDimension; d++)
    {
    sprintf(splitName, "%s%d%s", prefix, d, dotPtr);

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput(imageVector[d]);
    writer->SetFileName( splitName );
    try
      {
      writer->Update();
      }
    catch ( itk::ExceptionObject & e )
      {
      std::cerr << "Exception detected while generating displacement field" << splitName;
      std::cerr << " : "  << e.GetDescription();
      }
    }

}

template< class TVector>
void
PicslImageHelper
::CopyWithMax(TVector &maxVec, const TVector &newVec)
{
for (unsigned int d=0; d<TVector::Dimension; d++)
    {
    if (maxVec[d] < newVec[d])
      {
      maxVec[d] = newVec[d];
      }
    }
}

template< class TVector>
void
PicslImageHelper
::CopyWithMin(TVector &minVec, const TVector &newVec)
{
  for (unsigned int d=0; d<TVector::Dimension; d++)
    {
    if (minVec[d] > newVec[d])
      {
      minVec[d] = newVec[d];
      }
    }
}

}

#endif
