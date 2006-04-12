#ifndef __itkvanHerkGilWermanBaseImageFilter_h
#define __itkvanHerkGilWermanBaseImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkSize.h"
#include <vector>

namespace itk {

/** \class vanHerkGilWermanBaseImageFilter 
 * \brief Implements the van Herk/Gil-Werman algorithm for erosions
 * and dilations on flat, linear structuring elements. This algorithm
 * has complexity independent of structuring element size.
 *
 * \author Richard Beare. Department of Medicine, Monash University,
 * Melbourne, Australia.
 *
 * \sa 
 * \ingroup MathematicalMorphologyImageFilters
 */

template<class TInputImage, class TOutputImage, 
	 class TFunction1>
class ITK_EXPORT vanHerkGilWermanBaseImageFilter : 
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef vanHerkGilWermanBaseImageFilter Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>
  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage InputImageType;
  typedef TOutputImage OutputImageType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::RegionType      InputImageRegionType;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef typename InputImageType::SizeType        ISizeType;
  typedef typename InputImageType::OffsetType      IOffsetType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;
  typedef typename OutputImageType::RegionType     OutputImageRegionType;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;
  typedef typename OutputImageType::SizeType       OSizeType;
  
  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(vanHerkGilWermanBaseImageFilter, 
               ImageToImageFilter);

  /**
   * Set/Get the value used for boundary conditions
   */
  itkSetMacro(BoundaryValue, typename TInputImage::PixelType);
  itkGetConstReferenceMacro(BoundaryValue, typename TInputImage::PixelType);

  /**
   *  Set/Get for kernel size
   */
  itkSetMacro(KernelLength, int);
  itkGetConstReferenceMacro(KernelLength, int);

  /**
   *  Set/Get for kernel orientation
   */
  itkSetMacro(Direction, int);
  itkGetConstReferenceMacro(Direction, int);

protected:
  vanHerkGilWermanBaseImageFilter();
  ~vanHerkGilWermanBaseImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** vanHerkGilWermanBaseImageFilter needs the entire input be
   * available. Thus, it needs to provide an implementation of
   * GenerateInputRequestedRegion(). */
  void GenerateInputRequestedRegion() ;

  /** vanHerkGilWermanBaseImageFilter will produce the entire output. */
  void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));
  
  void GenerateData();
  

private:
  vanHerkGilWermanBaseImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  typename TInputImage::PixelType m_BoundaryValue;
  // the kernel size - must be odd
  int m_KernelLength;
  int m_Direction;

  typedef typename OutputImageType::IndexType OutIndexType;
  typedef typename InputImageType::IndexType InIndexType;

  // Line iterator
  typedef itk::ImageLinearConstIteratorWithIndex<InputImageType> InputLineIteratorType;
  // type for the running extrema
  typedef typename std::vector<InputImagePixelType> RunningExtType;

#if 0
  void forwardMax(InputLineIteratorType LineIt, RunningExtType &runningBuf);
  void backwardMax(InputLineIteratorType LineIt, RunningExtType &runningBuf);

  void forwardMaxBC(InputLineIteratorType LineIt, int middle, 
		    int LineLength, RunningExtType &runningBuf);
  void backwardMaxBC(InputLineIteratorType LineIt, RunningExtType &runningBuf);
#else

  typedef typename InputImageType::OffsetType InOffsetType;
  typedef typename std::vector<InOffsetType> InOffsetVecType;

  void mkOffsets(const int direction, InOffsetVecType &result);
  void forwardExtreme(const InIndexType start, const InOffsetVecType offsets, 
		      typename TInputImage::ConstPointer input,
		      const int blocks,
		      RunningExtType &runningBuf, RunningExtType &lineBuf);

  void reverseExtreme(const RunningExtType lineBuf, //const int start_last_block,
		      const int blocks, RunningExtType &runningBuf);


#endif
  TFunction1 Extreme;
} ; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkvanHerkGilWermanBaseImageFilter.txx"
#endif


#endif
