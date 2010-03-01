/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTConvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTConvolutionImageFilter_txx
#define __itkFFTConvolutionImageFilter_txx

#include "itkFFTConvolutionImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkFFTZeroPaddingImageFilter.h"
#include "itkNormalizeToConstantImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.h"
#include "itkRegionFromReferenceImageFilter.h"

namespace itk {

template <class TInputImage, class TKernelImage, class TOutputImage, class TFFTPrecision>
FFTConvolutionImageFilter<TInputImage, TKernelImage, TOutputImage, TFFTPrecision>
::FFTConvolutionImageFilter()
{
  this->SetNumberOfRequiredInputs(2);
}

template <class TInputImage, class TKernelImage, class TOutputImage, class TFFTPrecision>
void 
FFTConvolutionImageFilter<TInputImage, TKernelImage, TOutputImage, TFFTPrecision>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  InputImageType * input0 = const_cast<InputImageType *>(this->GetInput());
  InputImageType * input1 = const_cast<KernelImageType *>(this->GetKernelImage());
  if ( !input0 || !input1 )
    { 
    return;
    }
  input0->SetRequestedRegion( input0->GetLargestPossibleRegion() );
  input1->SetRequestedRegion( input1->GetLargestPossibleRegion() );
}


template<class TInputImage, class TKernelImage, class TOutputImage, class TFFTPrecision>
void
FFTConvolutionImageFilter<TInputImage, TKernelImage, TOutputImage, TFFTPrecision>
::GenerateData()
{
  this->AllocateOutputs();
  const InputImageType * input = this->GetInput();
  const KernelImageType * kernel = this->GetKernelImage();
  OutputImageType * output = this->GetOutput();

  typedef typename itk::Image< FFTPrecisionType, ImageDimension > InternalImageType;
  
  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  typedef itk::NormalizeToConstantImageFilter< KernelImageType, InternalImageType > NormType;
  typename NormType::Pointer norm = NormType::New();
  norm->SetInput( kernel );
  norm->SetNumberOfThreads( this->GetNumberOfThreads() );
  progress->RegisterInternalFilter( norm, 0.01f );

  typedef itk::FFTZeroPaddingImageFilter< InputImageType, InternalImageType, InternalImageType, InternalImageType > PadType;
  typename PadType::Pointer pad = PadType::New();
  pad->SetInput( input );
  pad->SetInputKernel( norm->GetOutput() );
  // TODO: check if this is needed
  pad->SetPadToPowerOfTwo( true );
  pad->SetNumberOfThreads( this->GetNumberOfThreads() );
  progress->RegisterInternalFilter( pad, 0.05f );

  typedef itk::FFTRealToComplexConjugateImageFilter< FFTPrecisionType, ImageDimension > FFTType;
  typename FFTType::Pointer fft = FFTType::New();
  fft->SetInput( pad->GetOutput() );
  fft->SetNumberOfThreads( this->GetNumberOfThreads() );
  progress->RegisterInternalFilter( fft, 0.25f );

  typedef itk::FFTShiftImageFilter< InternalImageType, InternalImageType > ShiftType;
  typename ShiftType::Pointer shift = ShiftType::New();
  shift->SetInput( pad->GetOutputKernel() );
  shift->SetNumberOfThreads( this->GetNumberOfThreads() );
  progress->RegisterInternalFilter( shift, 0.04f );
  
  typename FFTType::Pointer fftk = FFTType::New();
  fftk->SetInput( shift->GetOutput() );
  fftk->SetNumberOfThreads( this->GetNumberOfThreads() );
  progress->RegisterInternalFilter( fftk, 0.25f );
  
  typedef itk::MultiplyImageFilter< typename FFTType::OutputImageType,
                                    typename FFTType::OutputImageType,
                                    typename FFTType::OutputImageType > MultType;
  typename MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, fft->GetOutput() );
  mult->SetInput( 1, fftk->GetOutput() );
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  progress->RegisterInternalFilter( mult, 0.1f );
  
  typedef itk::FFTComplexConjugateToRealImageFilter< FFTPrecisionType, ImageDimension > IFFTType;
  typename IFFTType::Pointer ifft = IFFTType::New();
  ifft->SetInput( mult->GetOutput() );
  ifft->SetNumberOfThreads( this->GetNumberOfThreads() );
  progress->RegisterInternalFilter( ifft, 0.25f );
  
  typedef itk::RegionFromReferenceImageFilter< InternalImageType, OutputImageType > CropType;
  typename CropType::Pointer crop = CropType::New();
  crop->SetInput( ifft->GetOutput() );
  crop->SetReferenceImage( input );
  crop->SetNumberOfThreads( this->GetNumberOfThreads() );
  progress->RegisterInternalFilter( crop, 0.05f );
  
  crop->GraftOutput( output );
  crop->Update();
  this->GraftOutput( crop->GetOutput() );
}
  
}// end namespace itk
#endif