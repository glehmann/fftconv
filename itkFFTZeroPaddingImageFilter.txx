/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTZeroPaddingImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTZeroPaddingImageFilter_txx
#define __itkFFTZeroPaddingImageFilter_txx

#include "itkFFTZeroPaddingImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkNumericTraits.h"
#include "itkMath.h"
#include "itkConstantPadImageFilter.h"
#include "itkChangeInformationImageFilter.h"

namespace itk {

template <class TInputImage, class TOutputImage>
FFTZeroPaddingImageFilter<TInputImage, TOutputImage>
::FFTZeroPaddingImageFilter()
{
  m_PadToPowerOfTwo = false;
  this->SetNumberOfRequiredInputs(2);
  this->SetNumberOfRequiredOutputs(2);
  this->SetNthOutput( 1, OutputImageType::New() );
}

template <class TInputImage, class TOutputImage>
void 
FFTZeroPaddingImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  InputImageType * input0 = const_cast<InputImageType *>(this->GetInput(0));
  InputImageType * input1 = const_cast<InputImageType *>(this->GetInput(1));
  if ( !input0 || !input1 )
    { 
    return;
    }
  
  OutputImageType * output = this->GetOutput();
  
  RegionType region = output->GetRequestedRegion();
  region.Crop( input0->GetLargestPossibleRegion() );
  input0->SetRequestedRegion( region );
  
  region = output->GetRequestedRegion();
  region.Crop( input1->GetLargestPossibleRegion() );
  input1->SetRequestedRegion( region );
}


template <class TInputImage, class TOutputImage>
void 
FFTZeroPaddingImageFilter<TInputImage, TOutputImage>
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();
  
  const InputImageType * input0 = this->GetInput(0);
  const InputImageType * input1 = this->GetInput(1);
  if ( !input0 || !input1 )
    { 
    return;
    }
  
  OutputImageType * output0 = this->GetOutput(0);
  OutputImageType * output1 = this->GetOutput(1);
  
  RegionType region0 = input0->GetLargestPossibleRegion();
  RegionType region1 = input1->GetLargestPossibleRegion();
  
  // increase the size of the output by twice the size of the kernel (there are 2 borders on each dim).
  SizeType size;
  IndexType idx;
  for( int i=0; i<ImageDimension; i++ )
    {
    size[i] = region0.GetSize()[i] + 2 * region1.GetSize()[i];
    idx[i] = region0.GetIndex()[i] - region1.GetSize()[i];
    if( m_PadToPowerOfTwo )
      {
      // we may have to change the size of the region so that it is a power of 2
      // lets find the closest power of two which is greater or equal to the actual size
      unsigned long s2 = (unsigned long)vcl_pow(2.0f, (float)Math::Ceil(vcl_log((float)size[i])/vcl_log(2.0)));
      idx[i] -= ( s2 - size[i] ) / 2;
      size[i] = s2;
      }
    }
  RegionType region( idx, size );
  output0->SetLargestPossibleRegion( region );
  output1->SetLargestPossibleRegion( region );
  
}


template<class TInputImage, class TOutputImage>
void
FFTZeroPaddingImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  this->AllocateOutputs();
  const InputImageType * input0 = this->GetInput(0);
  const InputImageType * input1 = this->GetInput(1);
  OutputImageType * output0 = this->GetOutput(0);
  OutputImageType * output1 = this->GetOutput(1);
  RegionType ir0 = input0->GetLargestPossibleRegion();
  RegionType ir1 = input1->GetLargestPossibleRegion();
  RegionType or0 = output0->GetLargestPossibleRegion();
  RegionType or1 = output1->GetLargestPossibleRegion();

  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  typedef typename itk::ConstantPadImageFilter< InputImageType, OutputImageType > PadType;
  SizeType s;
  
  typename PadType::Pointer pad0 = PadType::New();
  pad0->SetInput( input0 );
  pad0->SetNumberOfThreads( this->GetNumberOfThreads() );
  for( int i=0; i<ImageDimension; i++ )
    {
    s[i] = ir0.GetIndex()[1] - or0.GetIndex()[i];
    }
  pad0->SetPadUpperBound( s );
  for( int i=0; i<ImageDimension; i++ )
    {
    s[i] = or0.GetSize()[i] - ( ir0.GetIndex()[1] - or0.GetIndex()[i] + ir0.GetSize()[i]);
    }
  pad0->SetPadLowerBound( s );
  progress->RegisterInternalFilter( pad0, 0.5f );
  pad0->GraftOutput( output0 );
  pad0->Update();
  this->GraftOutput( pad0->GetOutput() );

  typename PadType::Pointer pad1 = PadType::New();
  pad1->SetInput( input1 );
  pad1->SetNumberOfThreads( this->GetNumberOfThreads() );
  for( int i=0; i<ImageDimension; i++ )
    {
    s[i] = ( or1.GetSize()[i] - ir1.GetSize()[i] ) / 2;
    }
  pad1->SetPadUpperBound( s );
  for( int i=0; i<ImageDimension; i++ )
    {
    s[i] = itk::Math::Ceil(( or1.GetSize()[i] - ir1.GetSize()[i] ) / 2.0 );
    }
  pad1->SetPadLowerBound( s );
  progress->RegisterInternalFilter( pad1, 0.5f );
  
  typedef typename itk::ChangeInformationImageFilter< OutputImageType > ChangeType;
  typename ChangeType::Pointer change = ChangeType::New();
  change->SetInput( pad1->GetOutput() );
  change->SetUseReferenceImage( true );
  change->SetReferenceImage( output1 );
  change->SetChangeRegion( true );
  // no progress for change - it does almost nothing
  
  change->GraftOutput( output1 );
  change->Update();
  this->GraftNthOutput( 1, change->GetOutput() );
}


template<class TInputImage, class TOutputImage>
void
FFTZeroPaddingImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "PadToPowerOfTwo: "  << m_PadToPowerOfTwo << std::endl;
}
  
}// end namespace itk
#endif