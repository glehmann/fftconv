#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"

#include "itkFFTWRealToComplexConjugateImageFilter.h"
#include "itkFFTWComplexConjugateToRealImageFilter.h"
#include "itkFFTConvolutionImageFilter.h"
#include "itkTestingMacros.h"


int main(int argc, char * argv[])
{

  if( argc < 4 )
    {
    std::cerr << "usage: " << argv[0] << " intput kernel output [gpf]" << std::endl;
    std::cerr << " input: the input image" << std::endl;
    std::cerr << " output: the output image" << std::endl;
    std::cerr << " gpf: greatest prime factor of the size of the padded image" << std::endl;
    // std::cerr << "  : " << std::endl;
    exit(1);
    }

  const int dim = 2;
  
  typedef unsigned char PType;
  typedef itk::Image< PType, dim > IType;

  typedef itk::ImageFileReader< IType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  // itk::SimpleFilterWatcher watcher_reader(reader, "reader");

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  // itk::SimpleFilterWatcher watcher_reader2(reader2, "reader2");

  typedef itk::FFTConvolutionImageFilter< IType > FFTConvolutionType;
  FFTConvolutionType::Pointer conv = FFTConvolutionType::New();
  conv->SetInput( reader->GetOutput() );
  conv->SetKernelImage( reader2->GetOutput() );
  // test default value
  TEST_SET_GET_VALUE( 13, conv->GetGreatestPrimeFactor() );
  if( argc >= 5 )
    {
    conv->SetGreatestPrimeFactor( atoi(argv[4]) );
    TEST_SET_GET_VALUE( atoi(argv[4]), conv->GetGreatestPrimeFactor() );
    }
  itk::SimpleFilterWatcher watcher_norm(conv, "conv");

  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( conv->GetOutput() );
  writer->SetFileName( argv[3] );
  // itk::SimpleFilterWatcher watcher_writer(writer, "writer");
  
  writer->Update();

  return 0;
}
