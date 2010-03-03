#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"

#include "itkFFTZeroPaddingImageFilter.h"
#include "itkTestingMacros.h"


int main(int argc, char * argv[])
{

  if( argc != 6 )
    {
    std::cerr << "usage: " << argv[0] << " intput inputkernel output outputkernel gpf" << std::endl;
    std::cerr << " input: the input image" << std::endl;
    std::cerr << " inputkernel: the input kernel image" << std::endl;
    std::cerr << " output: the padded output image" << std::endl;
    std::cerr << " outputkernel: the padded output kernel image" << std::endl;
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

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );

  typedef itk::FFTZeroPaddingImageFilter< IType, IType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetInput( 1, reader2->GetOutput() );
  
  // test default value
  TEST_SET_GET_VALUE( 13, filter->GetGreatestPrimeFactor() );
  
  // test PadToPowerOfTwo behavior
  TEST_SET_GET_VALUE( 0, filter->GetPadToPowerOfTwo() );
  filter->SetPadToPowerOfTwo( true );
  TEST_SET_GET_VALUE( 1, filter->GetPadToPowerOfTwo() );
  TEST_SET_GET_VALUE( 2, filter->GetGreatestPrimeFactor() );
  
  // ok, now use the user provided value
  filter->SetGreatestPrimeFactor( atoi(argv[5]) );
  TEST_SET_GET_VALUE( atoi(argv[5]), filter->GetGreatestPrimeFactor() );

  itk::SimpleFilterWatcher watcher(filter, "filter");

  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();
  
  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetInput( filter->GetOutput(1) );
  writer2->SetFileName( argv[4] );
  writer2->Update();

  return 0;
}

