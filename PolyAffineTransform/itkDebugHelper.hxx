#ifndef __itkDebugHelper_hxx
#define __itkDebugHelper_hxx

#include "itkDebugHelper.h"

namespace itk
{

DebugHelper
::DebugHelper()
{
}

DebugHelper
::~DebugHelper()
{
}

void
DebugHelper::PrintSelf(std::ostream & os, Indent indent) const
{
}

template< class TImage>
static void 
DebugHelper
::WriteImage(typename TImage::Pointer image, char *fname)
{
  if (!DebugHelper::m_Debug)
    {
    return;
    }
  typedef itk::ImageFileWriter<TImage> WriterType;
  
  WriterType::Pointer writer = WriterType::New();
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
static void 
DebugHelper
::WriteDisplacementField(typename TField::Pointer field, char *fname)
{
  if (!DebugHelper::m_Debug)
    {
    return;
    }
  const unsigned int ImageDimension = TField::ImageDimension;
  typedef TField::PixelType VectorType;
  typedef TField::IndexType IndexType;
  typedef VectorType::ValueType ScalarType;  

  typedef itk::Image<ScalarType, ImageDimension> ScalarImageType;
  typedef itk::ImageFileWriter<ScalarImageType> WriterType;

  char *dotPtr = strrchr(fname, '.');
  int  dotPos = dotPtr - fname;
  char prefix[256], splitName[256];

  strcpy_s(prefix, fname);
  prefix[dotPos] = '\0';

  typedef ImageRegionIteratorWithIndex< TField > IteratorType;
  IteratorType it( field, field->GetLargestPossibleRegion() );

  VectorType vec;
  IndexType index;
  
  std::vector<ScalarImageType::Pointer> imageVector(ImageDimension);
  for (unsigned int i=0; i<ImageDimension; i++)
    {
    ScalarImageType::Pointer image = ScalarImageType::New();
    imageVector[i] = image;
    imageVector[i]->SetRegions(field->GetLargestPossibleRegion());
    imageVector[i]->CopyInformation(field);
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
    sprintf_s(splitName, "%s%d%s", prefix, d, dotPtr);

    WriterType::Pointer writer = WriterType::New();
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
}

#endif
