/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTZeroPaddingImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-01-28 18:14:36 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTZeroPaddingImageFilter_h
#define __itkFFTZeroPaddingImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkConceptChecking.h"

namespace itk {

/** \class FFTZeroPaddingImageFilter
 * \brief Produce a binary image where foreground is the regional maxima of the
 * input image
 *
 * Regional maxima are flat zones surounded by pixels of lower value.
 *
 * If the input image is constant, the entire image can be considered as a
 * maxima or not.  The desired behavior can be selected with the
 * SetFlatIsMaxima() method.
 * 
 * \author Gaetan Lehmann
 *
 * This class was contributed to the Insight Journal by author Gaetan Lehmann.
 * Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas,
 * France. The paper can be found at
 * http://insight-journal.org/midas/handle.php?handle=1926/153
 *
 * \sa ValuedFFTZeroPaddingImageFilter
 * \sa HConvexImageFilter 
 * \sa RegionalMinimaImageFilter
 *
 * \ingroup MathematicalMorphologyImageFilters
 */
template<class TInputImage, class TOutputImage>
class ITK_EXPORT FFTZeroPaddingImageFilter : 
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FFTZeroPaddingImageFilter Self;

  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;

  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef TOutputImage                             OutputImageType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;
  typedef typename InputImageType::RegionType      RegionType;
  typedef typename InputImageType::IndexType       IndexType;
  typedef typename InputImageType::SizeType        SizeType;
  
  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(FFTZeroPaddingImageFilter, ImageToImageFilter);

  /**
   * Set/Get whether the images must be padded to a size equal to a power of two.
   * This is required for vnl implementation of FFT, but not for FFTW.
   * The default is false.
   */
  itkSetMacro(PadToPowerOfTwo, bool);
  itkGetConstMacro(PadToPowerOfTwo, bool);
  itkBooleanMacro(PadToPowerOfTwo);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasPixelTraitsCheck,
    (Concept::HasPixelTraits<InputImagePixelType>));
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImagePixelType>));
  /** End concept checking */
#endif


protected:
  FFTZeroPaddingImageFilter();
  ~FFTZeroPaddingImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateInputRequestedRegion();
  void GenerateOutputInformation();
  
  /** Single-threaded version of GenerateData.  This filter delegates
   * to other filters. */
  void GenerateData();
  

private:
  FFTZeroPaddingImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool m_PadToPowerOfTwo;

}; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFFTZeroPaddingImageFilter.txx"
#endif

#endif
