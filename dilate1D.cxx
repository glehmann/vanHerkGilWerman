#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBasicDilateImageFilter.h"
#include "itkvanHerkGilWermanBaseImageFilter.h"
#include "itkNeighborhood.h"
#include "itkTimeProbe.h"
#include <vector>
#include "itkCommand.h"
#include "itkNumericTraits.h"

template < class TFilter >
class ProgressCallback : public itk::Command
{
public:
  typedef ProgressCallback   Self;
  typedef itk::Command  Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  itkTypeMacro( IterationCallback, Superclass );
  itkNewMacro( Self );

  /** Type defining the optimizer. */
  typedef    TFilter     FilterType;

  /** Method to specify the optimizer. */
  void SetFilter( FilterType * filter )
    {
    m_Filter = filter;
    m_Filter->AddObserver( itk::ProgressEvent(), this );
    }

  /** Execute method will print data at each iteration */
  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object *, const itk::EventObject & event)
    {
    std::cout << m_Filter->GetNameOfClass() << ": " << m_Filter->GetProgress() << std::endl;
    }

protected:
  ProgressCallback() {};
  itk::WeakPointer<FilterType>   m_Filter;
};

// a max functor, until I can figure out how to use the one in std
template <class pixtype>
class MaxFunctor
{
public:
  MaxFunctor(){}
  ~MaxFunctor(){}
  inline pixtype operator()(const pixtype &A, const pixtype &B)
  {
    return std::max(A, B);
  }
};

int main(int, char * argv[])
{
  const int dim = 2;
  typedef unsigned char PType;
  typedef itk::Image< PType, dim >    IType;
  
  // read the input image
  typedef itk::ImageFileReader< IType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New(); 
  reader->SetFileName( argv[1] );

  int radius = atoi(argv[4]);

#if 1
  typedef itk::Neighborhood<PType, dim> SRType;
  SRType kernel;
  SRType::SizeType ksize;
  ksize[0]=radius;
  ksize[1]=0;
  kernel.SetRadius(ksize);
  for( SRType::Iterator kit=kernel.Begin(); kit!=kernel.End(); kit++ )
    {
    *kit = 1;
    }

  typedef itk::BasicDilateImageFilter< IType, IType, SRType > BDilateType;
  BDilateType::Pointer dilate = BDilateType::New();
  dilate->SetInput( reader->GetOutput() );
  dilate->SetKernel( kernel );
#endif

  typedef itk::vanHerkGilWermanBaseImageFilter<IType, IType, MaxFunctor<PType> > DilateType;

  DilateType::Pointer vDilate = DilateType::New();
  vDilate->SetInput(reader->GetOutput());
  vDilate->SetKernelLength(2*radius+1);
  vDilate->SetBoundaryValue(itk::NumericTraits<PType>::NonpositiveMin());
  // write 
  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer vHwriter = WriterType::New();
  vHwriter->SetInput( vDilate->GetOutput() );
  vHwriter->SetFileName( argv[2] );
  vHwriter->Update();
#if 1  
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( dilate->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();
#endif
  
  return 0;
}

