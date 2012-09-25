#ifndef __itkDebugHelper_h
#define __itkDebugHelper_h

#include <string.h>

#include "itkMacro.h"
#include "itkObject.h"
#include "itkSmartPointer.h"

namespace itk
{

class DebugHelper:
  public Object
{
public:
  /** Standard typedefs   */
  typedef DebugHelper                   Self;
  typedef Object                        Superclass;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(DebugHelper, Object);
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  template< class TField>
  static void WriteDisplacementField(typename TField::Pointer field, char *fname);

  template< class TImage>
  static void WriteImage(typename TImage::Pointer image, char *fname);

  static const bool m_Debug = true;

protected:
  DebugHelper();
  virtual ~DebugHelper();

  /** Print contents of an DebugHelper */
  void PrintSelf(std::ostream & s, Indent indent) const;

private:

  DebugHelper(const Self & other);
  const Self & operator=(const Self &);

};

}

#if ITK_TEMPLATE_TXX
#include "itkDebugHelper.hxx"
#endif

#endif /* __itkDebugHelper_h */
