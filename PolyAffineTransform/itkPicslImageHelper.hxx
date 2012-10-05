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

bool PicslImageHelper::m_Debug = true;
char PicslImageHelper::m_FilePath[256] = {'\0'};

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

char * 
PicslImageHelper
::AppendPathToFileName(char *fullFileName, const char *fname)
{
  // USED TO FILE SYSTEM OPERATION
  struct stat st;
  char pathSeparator = '/';
  if(stat(m_FilePath, &st) != 0)
    {
    #if ( defined( _WIN64 ) || defined( _WIN32 ) || defined( WIN32 ) )
    pathSeparator = '\\';
    if (!CreateDirectory(m_FilePath, NULL))
      {
            std::string msg("The output directory did not exists, can't be created\n");
            std::cerr << msg;
            //throw std::runtime_error(msg);
      }
    #else
    pathSeparator = '/';
    if(mkdir(m_FilePath, S_IRWXU) != 0)
      {
            std::string msg("The output directory did not exist, can't be created\n");
            std::cerr << msg;
            //throw std::runtime_error(msg);
      }
    #endif
    }
  if (strlen(m_FilePath) == 0)
    {
    sprintf(fullFileName, "%s%s", m_FilePath, fname);
    }
  else
    {
    sprintf(fullFileName, "%s%c%s", m_FilePath, pathSeparator, fname);
    }
  return fullFileName;
}

char * 
PicslImageHelper
::AppendNumberToFileName(char *numberedFileName, const char *fname, int id)
{
  const char *dotPtr = strchr(fname, '.');
  int  dotPos = dotPtr - fname;
  char prefix[256], splitName[256];

  strcpy(prefix, fname);
  prefix[dotPos] = '\0';

  sprintf(splitName, "%s%d%s", prefix, id, dotPtr);
  strcpy(numberedFileName, splitName);

  return numberedFileName;
}

template< class TImage>
void 
PicslImageHelper
::WriteImage(typename TImage::Pointer image, char *fname, int id)
{
  char newName[256];
  PicslImageHelper::AppendNumberToFileName(newName, fname, id);
  PicslImageHelper::WriteImage<TImage>(image, newName);
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
  char fullName[256];
  PicslImageHelper::AppendPathToFileName(fullName, fname);

  typedef ImageFileWriter<TImage> WriterType;
  
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(image);
  writer->SetFileName( fullName );

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
::WriteDisplacementField(typename TField::Pointer field, char *fname, int id)
{
  char newName[256];
  PicslImageHelper::AppendNumberToFileName(newName, fname, id);
  PicslImageHelper::WriteDisplacementField<TField>(field, newName);
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
    char newName[256];
    PicslImageHelper::AppendNumberToFileName(newName, fname, d);
    char fullName[256];
    PicslImageHelper::AppendPathToFileName(fullName, newName);

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput(imageVector[d]);
    writer->SetFileName( fullName );
    try
      {
      writer->Update();
      }
    catch ( itk::ExceptionObject & e )
      {
      std::cerr << "Exception detected while generating displacement field" << newName;
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
