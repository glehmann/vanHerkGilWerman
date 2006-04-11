#ifndef __itkvanHerkGilWermanBaseImageFilter_txx
#define __itkvanHerkGilWermanBaseImageFilter_txx

#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"
#include "itkvanHerkGilWermanBaseImageFilter.h"
#include "itkProgressReporter.h"
#include "itkImageLinearIteratorWithIndex.h"
//#include "itkImageRegionIterator.h"
#include <iomanip>

namespace itk {

template <class TInputImage, class TOutputImage, class TFunction1>
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::vanHerkGilWermanBaseImageFilter()
{
  m_KernelLength=3;
  m_Direction = 0;
}

template <class TInputImage, class TOutputImage, class TFunction1>
void 
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // We need all the input.
  InputImagePointer input = const_cast<InputImageType *>(this->GetInput());
  if( !input )
    { return; }

  input->SetRequestedRegion( input->GetLargestPossibleRegion() );
}


template <class TInputImage, class TOutputImage, class TFunction1>
void 
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::EnlargeOutputRequestedRegion(DataObject *)
{
  this->GetOutput()
    ->SetRequestedRegion( this->GetOutput()->GetLargestPossibleRegion() );
}

#if 0
template<class TInputImage, class TOutputImage, class TFunction1>
void
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::GenerateData()
{
  // This version implements the vanHerk algorithm on a per block
  // basis, which means that it can't be used in place, but it doesn't
  // require minimal buffer and should have good cache characteristics.
  // Allocate the output
  this->AllocateOutputs();
  typename TOutputImage::Pointer output = this->GetOutput();
  typename TInputImage::ConstPointer input = this->GetInput();
  // the buffers
  RunningExtType ForwardRunningExt(m_KernelLength);
  RunningExtType BackwardRunningExt(m_KernelLength);

  int direction = 0;

  typedef itk::ImageLinearIteratorWithIndex<OutputImageType> OutputLineIteratorType;
  OutputLineIteratorType outLineIt(output, output->GetRequestedRegion());
  InputLineIteratorType inLineIt(input, output->GetRequestedRegion());
  OSizeType OSize = output->GetRequestedRegion().GetSize();
  int LineLength  = OSize[direction];
  inLineIt.SetDirection(direction);
  outLineIt.SetDirection(direction);

  int hKernel    = m_KernelLength/2;
  int halfBlocks = LineLength/hKernel;
  IOffsetType MOffset;
  MOffset.Fill(0);
  output->FillBuffer(0);
  for (inLineIt.GoToBegin(), outLineIt.GoToBegin(); 
       ! inLineIt.IsAtEnd(); 
       inLineIt.NextLine(), outLineIt.NextLine())
    {
    inLineIt.GoToBeginOfLine();
    outLineIt.GoToBeginOfLine();
    // the algorithm works by creating forward and backward running
    // extrema from points every KernelLength along the line of
    // interest, starting at KernelLength/2
    InputLineIteratorType middleIt = inLineIt;
    // do the first few pixels for which the backward running extrema
    // will run over the boundary
    forwardMax(middleIt, ForwardRunningExt);
    backwardMaxBC(middleIt, BackwardRunningExt);
    for (int k=hKernel;k<m_KernelLength-1;++k, ++outLineIt) 
      {
      //std::cout << k << std::endl;
      outLineIt.Set(Extreme(ForwardRunningExt[k], BackwardRunningExt[k]));
      }
    int middle;
    for (middle = m_KernelLength-1;
	 middle < LineLength - 3*hKernel;
	 middle += m_KernelLength)
      {
      MOffset[direction] = middle;
      middleIt.SetIndex(inLineIt.GetIndex() + MOffset);
      //std::cout << middle << " " << LineLength << std::endl;
      forwardMax(middleIt, ForwardRunningExt);
      backwardMax(middleIt, BackwardRunningExt);
      MOffset[direction] = middle - hKernel;
      outLineIt.SetIndex(inLineIt.GetIndex() + MOffset);
      for (int k=0;k<m_KernelLength;++k, ++outLineIt)
	{
	//std::cout << middle - hKernel + k << " " << (int)ForwardRunningExt[k] << " " << (int)BackwardRunningExt[k] << std::endl;
	outLineIt.Set(Extreme(ForwardRunningExt[k], BackwardRunningExt[k]));
	}
      }
    // Now for the last few pixels
    //std::cout << "Exit " << middle << " " << LineLength << std::endl;
    MOffset[direction] = middle;
    middleIt.SetIndex(inLineIt.GetIndex() + MOffset);
    forwardMaxBC(middleIt, middle, LineLength, ForwardRunningExt);
    backwardMax(middleIt, BackwardRunningExt);
    for (int k=0, i = middle - hKernel; i < LineLength; ++i, ++k, ++outLineIt)
      {
      outLineIt.Set(Extreme(ForwardRunningExt[k], BackwardRunningExt[k]));
      }
    //break;
    }
  

}

template<class TInputImage, class TOutputImage, class TFunction1>
void
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::forwardMax(InputLineIteratorType LineIt, RunningExtType &runningBuf)
{
  // Version that does not need a bounds check
  InputImagePixelType Ext = m_BoundaryValue;
  for (int i = 0;i < m_KernelLength; i++, ++LineIt)
    {
    Ext = Extreme(Ext, LineIt.Get());
    runningBuf[i] = Ext;
    }
}

template<class TInputImage, class TOutputImage, class TFunction1>
void
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::forwardMaxBC(InputLineIteratorType LineIt, int middle, 
	       int LineLength, RunningExtType &runningBuf)
{
  // Version that does not need a bounds check
  InputImagePixelType Ext = m_BoundaryValue;
  int j=0;
  for (int i = middle;i < LineLength; i++, ++LineIt, j++)
    {
    Ext = Extreme(Ext, LineIt.Get());
    runningBuf[j] = Ext;
    }
  // copy the last value to the end
  while (j < m_KernelLength)
    {
    runningBuf[j] = Ext;
    j++;
    }
}

template<class TInputImage, class TOutputImage, class TFunction1>
void
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::backwardMax(InputLineIteratorType LineIt, RunningExtType &runningBuf)
{
  // Version that deals with end of line bounds
  InputImagePixelType Ext = m_BoundaryValue;
  for (int i = m_KernelLength-1;i >= 0; i--, --LineIt)
    {
    Ext = Extreme(Ext, LineIt.Get());
    runningBuf[i] = Ext;
    }

}
template<class TInputImage, class TOutputImage, class TFunction1>
void
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::backwardMaxBC(InputLineIteratorType LineIt, RunningExtType &runningBuf)
{
  // Version that deals with bounds - only the first value is 
  // in the image, so just copy it
  InputImagePixelType Ext = LineIt.Get();
  for (int i = m_KernelLength-1;i >= 0; i--, --LineIt)
    {
    runningBuf[i] = Ext;
    }

}
#else

template<class TInputImage, class TOutputImage, class TFunction1>
void
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::GenerateData()
{
  this->AllocateOutputs();
  typename TOutputImage::Pointer output = this->GetOutput();
  typename TInputImage::ConstPointer input = this->GetInput();
  
  // create a vector of offsets defining the lines along which we will
  // operate
  InOffsetVecType Offsets;
  int direction = m_Direction;
  mkOffsets(direction, Offsets);

  // we want to make the size of the running extrema buffers
  // correspond to an integral number of blocks because it simplifies
  // the management
  int blocks = Offsets.size()/m_KernelLength;
  int leftovers = Offsets.size() % m_KernelLength;
  int runningBufSize, start_last_block;
  if (leftovers == 0) 
    {
    runningBufSize = Offsets.size();
    start_last_block = (blocks - 1) * m_KernelLength;
    }
  else
    {
    runningBufSize = m_KernelLength * (blocks + 1);
    start_last_block = blocks * m_KernelLength;
    }
//  std::cout << "leftovers=" << leftovers << " start_block = " << start_last_block << std::endl;
  // the buffers
  RunningExtType ForwardRunningExt(runningBufSize);
  RunningExtType BackwardRunningExt(runningBufSize);
  // this one will store the pixels so we don't need to get them from
  // the image a second time
  RunningExtType LineBuffer(Offsets.size()); 

  // use a line iterator to find the index at the beginning of each
  // line. Don't iterate along the line - we have already stored the
  // offsets. I'm not sure this is the most efficient way to go, but
  // it should be more readily adaptable to directions that aren't
  // parallel to the axes
  InputLineIteratorType inLineIt(this->GetInput(), this->GetOutput()->GetRequestedRegion());
  inLineIt.SetDirection(direction);
  inLineIt.GoToBegin();
  while (!inLineIt.IsAtEnd())
    {
    InIndexType start = inLineIt.GetIndex();

    forwardExtreme(start, Offsets, input, blocks, ForwardRunningExt, LineBuffer);
    reverseExtreme(LineBuffer, start_last_block, blocks, BackwardRunningExt);

    for (int i = 0, j=m_KernelLength/2;i < Offsets.size(); i++, j++)
      {
      output->SetPixel(start + Offsets[i], Extreme(ForwardRunningExt[j], BackwardRunningExt[j]));
      }

    inLineIt.NextLine();
    }
}

template<class TInputImage, class TOutputImage, class TFunction1>
void
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::mkOffsets(const int direction, InOffsetVecType &result)
{
  InputLineIteratorType inLineIt(this->GetInput(), this->GetOutput()->GetRequestedRegion());
  inLineIt.SetDirection(direction);
  inLineIt.GoToBegin();
  inLineIt.GoToBeginOfLine();
  InIndexType start = inLineIt.GetIndex();

  for (inLineIt.GoToBeginOfLine(); !inLineIt.IsAtEndOfLine(); ++inLineIt)
    {
    InOffsetType o = inLineIt.GetIndex() - start;
    result.push_back(o);
    }

//  return(result);
}

template<class TInputImage, class TOutputImage, class TFunction1>
void
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::forwardExtreme(const InIndexType start, const InOffsetVecType offsets,
		 typename TInputImage::ConstPointer input, 
		 const int blocks,
		 RunningExtType &runningBuf, RunningExtType &lineBuf)
{
  // generate the forward running max
  // also store the pixel values in the line buffer for quick access
  // when computing the reverse running max
  InputImagePixelType Ext = m_BoundaryValue;
#if 0
  for (int i = 0; i < offsets.size(); i++)
    {
    // restart the running extreme
    if ((i % m_KernelLength) == 0) 
      {
      Ext = m_BoundaryValue;
      }
    InputImagePixelType V = input->GetPixel(start + offsets[i]);
    lineBuf[i]=V;
    Ext = Extreme(Ext,V);
    runningBuf[i] = Ext;
    }
#else
  // slightly messier looking version that should minimize loop complexity
  int i = 0;
  for (int k=0;k<blocks;k++)
    {
    Ext = m_BoundaryValue;
    for (int j=0;j<m_KernelLength;j++, i++)
      {
      InputImagePixelType V = input->GetPixel(start + offsets[i]);
      lineBuf[i]=V;
      Ext = Extreme(Ext,V);
      runningBuf[i] = Ext;
      }
    }
  // finish the last of the block
  Ext = m_BoundaryValue;
  for (; i < offsets.size(); i++)
    {
    InputImagePixelType V = input->GetPixel(start + offsets[i]);
    lineBuf[i]=V;
    Ext = Extreme(Ext,V);
    runningBuf[i] = Ext;
    }
#endif
  // finish the rest
  for (int j = offsets.size();j < runningBuf.size();j++)
    {
    runningBuf[j] = Ext;
    }
}

template<class TInputImage, class TOutputImage, class TFunction1>
void
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::reverseExtreme(const RunningExtType lineBuf, const int start_last_block,
		 const int blocks,
		 RunningExtType &runningBuf)
{
  InputImagePixelType Ext = m_BoundaryValue;
  int LineIdx = start_last_block;
#if 0
  for (int i = runningBuf.size() - 1; i >= this->m_KernelLength - 1; i--, LineIdx--) 
    {
    if ((LineIdx ) % m_KernelLength == 0)
      {
      Ext = m_BoundaryValue;
      }
    Ext = Extreme(Ext, lineBuf[LineIdx]);
    runningBuf[i] = Ext;
    }
#else
  int i = runningBuf.size() - 1;
  for (int k=0;k<blocks;k++)
    {
    Ext = m_BoundaryValue;
    for (int j=0;j<m_KernelLength;j++, LineIdx--, i--)
      {
      Ext = Extreme(Ext, lineBuf[LineIdx]);
      runningBuf[i] = Ext;
      }
    }
  Ext = lineBuf[LineIdx];
#endif
  // now fill the start of the buffer -- this corresponds to the
  // boundary condition
  for (int j = 0; j < m_KernelLength-1;j++) 
    {
    runningBuf[j] = Ext;
    }
}

#endif
template<class TInputImage, class TOutputImage, class TFunction1>
void
vanHerkGilWermanBaseImageFilter<TInputImage, TOutputImage, TFunction1>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "KernelLength   : "  << m_KernelLength << std::endl;
  os << indent << "Boundary value : "  << m_BoundaryValue << std::endl;
}
  
} // end namespace itk

#endif
