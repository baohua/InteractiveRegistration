#ifndef __itkPicslImageHelper_h
#define __itkPicslImageHelper_h

#include <string.h>

#include "itkMacro.h"
#include "itkObject.h"
#include "itkSmartPointer.h"

namespace itk
{

class PicslImageHelper:
  public Object
{
public:
  /** Standard typedefs   */
  typedef PicslImageHelper                   Self;
  typedef Object                        Superclass;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(PicslImageHelper, Object);
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  template< class TField>
  static void WriteDisplacementField(typename TField::Pointer field, char *fname);

  template< class TImage>
  static void WriteImage(typename TImage::Pointer image, char *fname);

  template< class TRegion>
  static typename itk::VectorContainer<int, typename TRegion::IndexType>::Pointer 
    GetCorners( const TRegion &region);

  template< class TVector >
  static void CopyWithMax(TVector &maxVec, const TVector &newVec);
  template< class TVector >
  static void CopyWithMin(TVector &minVec, const TVector &newVec);

  static const bool m_Debug = true;

protected:
  PicslImageHelper();
  virtual ~PicslImageHelper();

  /** Print contents of an PicslImageHelper */
  void PrintSelf(std::ostream & s, Indent indent) const;

private:

  PicslImageHelper(const Self & other);
  const Self & operator=(const Self &);

};

}

#if ITK_TEMPLATE_TXX
#include "itkPicslImageHelper.hxx"
#endif

#endif /* __itkPicslImageHelper_h */
