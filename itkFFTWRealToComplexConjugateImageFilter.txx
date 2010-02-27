/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTWRealToComplexConjugateImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2010-02-26 23:50:55 $
  Version:   $Revision: 1.13 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTWRealToComplexConjugateImageFilter_txx
#define __itkFFTWRealToComplexConjugateImageFilter_txx

#include "itkFFTWRealToComplexConjugateImageFilter.h"
#include "itkFFTRealToComplexConjugateImageFilter.txx"
#include <iostream>
#include "itkIndent.h"
#include "itkMetaDataObject.h"
#include "itkProgressReporter.h"

namespace itk
{
/** TODO:  There should be compile time type checks so that
           if only USE_FFTWF is defined, then only floats are valid.
           and if USE_FFTWD is defined, then only doubles are valid.
*/

template <typename TPixel, unsigned int VDimension>
void
FFTWRealToComplexConjugateImageFilter<TPixel,VDimension>::
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

  const typename TInputImageType::SizeType&   inputSize
    = inputPtr->GetLargestPossibleRegion().GetSize();
  const typename TOutputImageType::SizeType&   outputSize
    = outputPtr->GetLargestPossibleRegion().GetSize();

  // figure out sizes
  // size of input and output aren't the same which is handled in the superclass,
  // sort of.
  // the input size and output size only differ in the fastest moving dimension
  unsigned int total_inputSize = 1;
  unsigned int total_outputSize = 1;

  for(unsigned i = 0; i < VDimension; i++)
    {
    total_inputSize *= inputSize[i];
    total_outputSize *= outputSize[i];
    }

    typename FFTWProxyType::PlanType plan;
    switch(VDimension)
      {
      case 1:
        plan = FFTWProxyType::Plan_dft_r2c_1d(inputSize[0],
                                             const_cast<TPixel*>(inputPtr->GetBufferPointer()),
                                             (typename FFTWProxyType::ComplexType*) outputPtr->GetBufferPointer(),
                                             FFTW_ESTIMATE);
        break;
      case 2:
        plan = FFTWProxyType::Plan_dft_r2c_2d(inputSize[1],
                                             inputSize[0],
                                             const_cast<TPixel*>(inputPtr->GetBufferPointer()),
                                             (typename FFTWProxyType::ComplexType*) outputPtr->GetBufferPointer(),
                                             FFTW_ESTIMATE);
        break;
      case 3:
        plan = FFTWProxyType::Plan_dft_r2c_3d(inputSize[2],
                                             inputSize[1],
                                             inputSize[0],
                                             const_cast<TPixel*>(inputPtr->GetBufferPointer()),
                                             (typename FFTWProxyType::ComplexType*)outputPtr->GetBufferPointer(),
//                                             FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
                                             FFTW_ESTIMATE);
        break;
      default:
        int *sizes = new int[VDimension];
        for(unsigned int i = 0; i < VDimension; i++)
          {
          sizes[(VDimension - 1) - i] = inputSize[i];
          }

        plan = FFTWProxyType::Plan_dft_r2c(VDimension,sizes,
                                          const_cast<TPixel*>(inputPtr->GetBufferPointer()),
                                          (typename FFTWProxyType::ComplexType*) outputPtr->GetBufferPointer(),
                                          FFTW_ESTIMATE);
        delete [] sizes;
      }
  FFTWProxyType::Execute(plan);
  FFTWProxyType::DestroyPlan(plan);
}

template <typename TPixel,unsigned int VDimension>
bool
FFTWRealToComplexConjugateImageFilter<TPixel,VDimension>::
FullMatrix()
{
  return false;
}

} // namespace itk

#endif //_itkFFTWRealToComplexConjugateImageFilter_txx
