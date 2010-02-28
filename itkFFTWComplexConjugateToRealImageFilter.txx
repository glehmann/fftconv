/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTWComplexConjugateToRealImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2010-02-26 23:50:55 $
  Version:   $Revision: 1.15 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTWComplexConjugateToRealImageFilter_txx
#define __itkFFTWComplexConjugateToRealImageFilter_txx

#include "itkFFTWComplexConjugateToRealImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.txx"
#include <iostream>
#include "itkIndent.h"
#include "itkMetaDataObject.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"

namespace itk
{

template <typename TPixel, unsigned int VDimension>
void
FFTWComplexConjugateToRealImageFilter<TPixel,VDimension>::
GenerateData()
{
  // get pointers to the input and output
  typename TInputImageType::ConstPointer  inputPtr  = this->GetInput();
  typename TOutputImageType::Pointer      outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // we don't have a nice progress to report, but at least this simple line
  // reports the begining and the end of the process
  ProgressReporter progress(this, 0, 1);

  // allocate output buffer memory
  outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
  outputPtr->Allocate();

  const typename TInputImageType::SizeType&   outputSize
    = outputPtr->GetLargestPossibleRegion().GetSize();
  const typename TOutputImageType::SizeType& inputSize
    = inputPtr->GetLargestPossibleRegion().GetSize();

  // figure out sizes
  // size of input and output aren't the same which is handled in the superclass,
  // sort of.
  // the input size and output size only differ in the fastest moving dimension
  unsigned int total_outputSize = 1;
  unsigned int total_inputSize = 1;

  for(unsigned i = 0; i < VDimension; i++)
    {
    total_outputSize *= outputSize[i];
    total_inputSize *= inputSize[i];
    }

    typename FFTWProxyType::ComplexType * in = new typename FFTWProxyType::ComplexType[total_inputSize];
    TPixel * out = outputPtr->GetBufferPointer();
    typename FFTWProxyType::PlanType plan;
    
    switch(VDimension)
      {
      case 1:
        plan = FFTWProxyType::Plan_dft_c2r_1d(outputSize[0],
                                       in,
				       out,
                                       FFTW_ESTIMATE,
				       this->GetNumberOfThreads());
        break;
      case 2:
        plan = FFTWProxyType::Plan_dft_c2r_2d(outputSize[1],outputSize[0],
                                       in,
				       out,
                                       FFTW_ESTIMATE,
				       this->GetNumberOfThreads());
        break;
      case 3:
        plan = FFTWProxyType::Plan_dft_c2r_3d(outputSize[2],outputSize[1],outputSize[0],
                                       in,
				       out,
                                       FFTW_ESTIMATE,
				       this->GetNumberOfThreads());
        break;
      default:
        int *sizes = new int[VDimension];
        for(unsigned int i = 0; i < VDimension; i++)
          {
          sizes[(VDimension - 1) - i] = outputSize[i];
          }
        plan = FFTWProxyType::Plan_dft_c2r(VDimension,sizes,
                                    in,
				    out,
                                    FFTW_ESTIMATE,
				    this->GetNumberOfThreads());
        delete [] sizes;
      }
  // copy the input, because it may be destroyed by computing the plan
  memcpy(in,
         inputPtr->GetBufferPointer(),
         total_inputSize * sizeof(typename FFTWProxyType::ComplexType));
  fftw::Proxy<TPixel>::Execute(plan);
  FFTWProxyType::DestroyPlan(plan);
  delete [] in;
  
  typedef ImageRegionIterator< TOutputImageType >   IteratorType;
  IteratorType it(outputPtr,outputPtr->GetLargestPossibleRegion());
  while( !it.IsAtEnd() )
    {
    it.Set( it.Value() / total_outputSize );
    ++it;
    }
}
template <typename TPixel,unsigned int VDimension>
bool
FFTWComplexConjugateToRealImageFilter<TPixel,VDimension>::
FullMatrix()
{
  return false;
}

}// namespace itk
#endif // _itkFFTWComplexConjugateToRealImageFilter_txx
