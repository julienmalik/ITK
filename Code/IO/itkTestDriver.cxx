/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
// define some itksys* things to make ShareForward.h happy
#define itksys_SHARED_FORWARD_DIR_BUILD ""
#define itksys_SHARED_FORWARD_PATH_BUILD ""
#define itksys_SHARED_FORWARD_PATH_INSTALL ""
#define itksys_SHARED_FORWARD_EXE_BUILD ""
#define itksys_SHARED_FORWARD_EXE_INSTALL ""

#include "itkWin32Header.h"
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "itkMultiThreader.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkDifferenceImageFilter.h"
#include "itksys/SystemTools.hxx"
// include SharedForward to avoid duplicating the code which find the library
// path variable
// name and the path separator
#include "itksys/SharedForward.h"
#include "itksys/Process.h"

#define ITK_TEST_DIMENSION_MAX 6

void usage()
{
  std::cerr << "usage: itkTestDriver [options] prg [args]" << std::endl;
  std::cerr << std::endl;
  std::cerr << "itkTestDriver alter the environment, run a test program and compare the images" << std::endl;
  std::cerr << "produced." << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "  --add-before-libpath PATH" << std::endl;
  std::cerr << "      Add a path to the library path environment. This option take care of" << std::endl;
  std::cerr << "      choosing the right environment variable for your system." << std::endl;
  std::cerr << "      This option can be used several times." << std::endl;
  std::cerr << std::endl;
  std::cerr << "  --add-before-env NAME VALUE" << std::endl;
  std::cerr << "      Add a VALUE to the variable name in the environment." << std::endl;
  std::cerr << "      The seperator used is the default one on the system." << std::endl;
  std::cerr << "      This option can be used several times." << std::endl;
  std::cerr << std::endl;
  std::cerr << "  --add-before-env-with-sep NAME VALUE SEP" << std::endl;
  std::cerr << "      Add a VALUE to the variable name in the environment using the provided separator." << std::endl;
  std::cerr << "      This option can be used several times." << std::endl;
  std::cerr << std::endl;
  std::cerr << "  --compare TEST BASELINE" << std::endl;
  std::cerr << "      Compare the TEST image to the BASELINE one." << std::endl;
  std::cerr << "      This option can be used several times." << std::endl;
  std::cerr << std::endl;
  std::cerr << "  --" << std::endl;
  std::cerr << "      The options after -- are not interpreted by this program and passed" << std::endl;
  std::cerr << "      directly to the test program." << std::endl;
  std::cerr << std::endl;
  std::cerr << "  --help" << std::endl;
  std::cerr << "      Display this message and exit." << std::endl;
  std::cerr << std::endl;
}

// Regression Testing Code

int RegressionTestImage(const char *testImageFilename, const char *baselineImageFilename, int reportErrors)
{
  // Use the factory mechanism to read the test and baseline files and convert
  // them to double
  typedef itk::Image< double, ITK_TEST_DIMENSION_MAX >        ImageType;
  typedef itk::Image< unsigned char, ITK_TEST_DIMENSION_MAX > OutputType;
  typedef itk::Image< unsigned char, 2 >                      DiffOutputType;
  typedef itk::ImageFileReader< ImageType >                   ReaderType;

  // Read the baseline file
  ReaderType::Pointer baselineReader = ReaderType::New();
  baselineReader->SetFileName(baselineImageFilename);
  try
    {
    baselineReader->UpdateLargestPossibleRegion();
    }
  catch ( itk::ExceptionObject & e )
    {
    std::cerr << "Exception detected while reading " << baselineImageFilename << " : "  << e.GetDescription();
    return 1000;
    }

  // Read the file generated by the test
  ReaderType::Pointer testReader = ReaderType::New();
  testReader->SetFileName(testImageFilename);
  try
    {
    testReader->UpdateLargestPossibleRegion();
    }
  catch ( itk::ExceptionObject & e )
    {
    std::cerr << "Exception detected while reading " << testImageFilename << " : "  << e.GetDescription() << std::endl;
    return 1000;
    }

  // The sizes of the baseline and test image must match
  ImageType::SizeType baselineSize;
  baselineSize = baselineReader->GetOutput()->GetLargestPossibleRegion().GetSize();
  ImageType::SizeType testSize;
  testSize = testReader->GetOutput()->GetLargestPossibleRegion().GetSize();

  if ( baselineSize != testSize )
    {
    std::cerr << "The size of the Baseline image and Test image do not match!" << std::endl;
    std::cerr << "Baseline image: " << baselineImageFilename
              << " has size " << baselineSize << std::endl;
    std::cerr << "Test image:     " << testImageFilename
              << " has size " << testSize << std::endl;
    return 1;
    }

  // Now compare the two images
  typedef itk::DifferenceImageFilter< ImageType, ImageType > DiffType;
  DiffType::Pointer diff = DiffType::New();
  diff->SetValidInput( baselineReader->GetOutput() );
  diff->SetTestInput( testReader->GetOutput() );
  diff->SetDifferenceThreshold(2.0);
  diff->UpdateLargestPossibleRegion();

  double status = diff->GetTotalDifference();

  // if there are discrepencies, create an diff image
  if ( status && reportErrors )
    {
    typedef itk::RescaleIntensityImageFilter< ImageType, OutputType > RescaleType;
    typedef itk::ExtractImageFilter< OutputType, DiffOutputType >     ExtractType;
    typedef itk::ImageFileWriter< DiffOutputType >                    WriterType;
    typedef itk::ImageRegion< ITK_TEST_DIMENSION_MAX >                RegionType;

    OutputType::IndexType index; index.Fill(0);
    OutputType::SizeType  size; size.Fill(0);

    RescaleType::Pointer rescale = RescaleType::New();
    rescale->SetOutputMinimum( itk::NumericTraits< unsigned char >::NonpositiveMin() );
    rescale->SetOutputMaximum( itk::NumericTraits< unsigned char >::max() );
    rescale->SetInput( diff->GetOutput() );
    rescale->UpdateLargestPossibleRegion();

    RegionType region;
    region.SetIndex(index);

    size = rescale->GetOutput()->GetLargestPossibleRegion().GetSize();
    for ( unsigned int i = 2; i < ITK_TEST_DIMENSION_MAX; i++ )
      {
      size[i] = 0;
      }
    region.SetSize(size);

    ExtractType::Pointer extract = ExtractType::New();
    extract->SetInput( rescale->GetOutput() );
    extract->SetExtractionRegion(region);

    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( extract->GetOutput() );

    std::cout << "<DartMeasurement name=\"ImageError\" type=\"numeric/double\">";
    std::cout << status;
    std::cout <<  "</DartMeasurement>" << std::endl;

    std::ostringstream diffName;
    diffName << testImageFilename << ".diff.png";
    try
      {
      rescale->SetInput( diff->GetOutput() );
      rescale->Update();
      }
    catch ( const std::exception & e )
      {
      std::cerr << "Error during rescale of " << diffName.str() << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch ( ... )
      {
      std::cerr << "Error during rescale of " << diffName.str() << std::endl;
      }
    writer->SetFileName( diffName.str().c_str() );
    try
      {
      writer->Update();
      }
    catch ( const std::exception & e )
      {
      std::cerr << "Error during write of " << diffName.str() << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch ( ... )
      {
      std::cerr << "Error during write of " << diffName.str() << std::endl;
      }

    std::cout << "<DartMeasurementFile name=\"DifferenceImage\" type=\"image/png\">";
    std::cout << diffName.str();
    std::cout << "</DartMeasurementFile>" << std::endl;

    std::ostringstream baseName;
    baseName << testImageFilename << ".base.png";
    try
      {
      rescale->SetInput( baselineReader->GetOutput() );
      rescale->Update();
      }
    catch ( const std::exception & e )
      {
      std::cerr << "Error during rescale of " << baseName.str() << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch ( ... )
      {
      std::cerr << "Error during rescale of " << baseName.str() << std::endl;
      }
    try
      {
      writer->SetFileName( baseName.str().c_str() );
      writer->Update();
      }
    catch ( const std::exception & e )
      {
      std::cerr << "Error during write of " << baseName.str() << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch ( ... )
      {
      std::cerr << "Error during write of " << baseName.str() << std::endl;
      }

    std::cout << "<DartMeasurementFile name=\"BaselineImage\" type=\"image/png\">";
    std::cout << baseName.str();
    std::cout << "</DartMeasurementFile>" << std::endl;

    std::ostringstream testName;
    testName << testImageFilename << ".test.png";
    try
      {
      rescale->SetInput( testReader->GetOutput() );
      rescale->Update();
      }
    catch ( const std::exception & e )
      {
      std::cerr << "Error during rescale of " << testName.str() << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch ( ... )
      {
      std::cerr << "Error during rescale of " << testName.str() << std::endl;
      }
    try
      {
      writer->SetFileName( testName.str().c_str() );
      writer->Update();
      }
    catch ( const std::exception & e )
      {
      std::cerr << "Error during write of " << testName.str() << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch ( ... )
      {
      std::cerr << "Error during write of " << testName.str() << std::endl;
      }

    std::cout << "<DartMeasurementFile name=\"TestImage\" type=\"image/png\">";
    std::cout << testName.str();
    std::cout << "</DartMeasurementFile>" << std::endl;
    }
  return ( status != 0 ) ? 1 : 0;
}

//
// Generate all of the possible baselines
// The possible baselines are generated fromn the baselineFilename using the
// following algorithm:
// 1) strip the suffix
// 2) append a digit .x
// 3) append the original suffix.
// It the file exists, increment x and continue
//
std::map< std::string, int > RegressionTestBaselines(char *baselineFilename)
{
  std::map< std::string, int > baselines;
  baselines[std::string(baselineFilename)] = 0;

  std::string originalBaseline(baselineFilename);

  int                    x = 0;
  std::string::size_type suffixPos = originalBaseline.rfind(".");
  std::string            suffix;
  if ( suffixPos != std::string::npos )
    {
    suffix = originalBaseline.substr( suffixPos, originalBaseline.length() );
    originalBaseline.erase( suffixPos, originalBaseline.length() );
    }
  while ( ++x )
    {
    std::ostringstream filename;
    filename << originalBaseline << "." << x << suffix;
    std::ifstream filestream( filename.str().c_str() );
    if ( !filestream )
      {
      break;
      }
    baselines[filename.str()] = 0;
    filestream.close();
    }
  return baselines;
}

int main(int ac, char *av[])
{
  std::vector< char * > args;
  typedef std::pair< char *, char * > ComparePairType;
  std::vector< ComparePairType > compareList;

  // parse the command line
  int  i = 1;
  bool skip = false;
  while ( i < ac )
    {
    if ( !skip && strcmp(av[i], "--add-before-libpath") == 0 )
      {
      if ( i + 1 >= ac )
        {
        usage();
        return 1;
        }
      std::string libpath = KWSYS_SHARED_FORWARD_LDPATH;
      libpath += "=";
      libpath += av[i + 1];
      char *oldenv = getenv(KWSYS_SHARED_FORWARD_LDPATH);
      if ( oldenv )
        {
        libpath += KWSYS_SHARED_FORWARD_PATH_SEP;
        libpath += oldenv;
        }
      itksys::SystemTools::PutEnv( libpath.c_str() );
      // on some 64 bit systems, LD_LIBRARY_PATH_64 is used before
      // LD_LIBRARY_PATH if it is set. It can lead the test to load
      // the system library instead of the expected one, so this
      // var must also be set
      if ( std::string(KWSYS_SHARED_FORWARD_LDPATH) == "LD_LIBRARY_PATH" )
        {
        std::string libpath = "LD_LIBRARY_PATH_64";
        libpath += "=";
        libpath += av[i + 1];
        char *oldenv = getenv("LD_LIBRARY_PATH_64");
        if ( oldenv )
          {
          libpath += KWSYS_SHARED_FORWARD_PATH_SEP;
          libpath += oldenv;
          }
        itksys::SystemTools::PutEnv( libpath.c_str() );
        }
      i += 2;
      }
    else if ( !skip && strcmp(av[i], "--add-before-env") == 0 )
      {
      if ( i + 2 >= ac )
        {
        usage();
        return 1;
        }
      std::string env = av[i + 1];
      env += "=";
      env += av[i + 2];
      char *oldenv = getenv(av[i + 1]);
      if ( oldenv )
        {
        env += KWSYS_SHARED_FORWARD_PATH_SEP;
        env += oldenv;
        }
      itksys::SystemTools::PutEnv( env.c_str() );
      i += 3;
      }
    else if ( !skip && strcmp(av[i], "--add-before-env-with-sep") == 0 )
      {
      if ( i + 3 >= ac )
        {
        usage();
        return 1;
        }
      std::string env = av[i + 1];
      env += "=";
      env += av[i + 2];
      char *oldenv = getenv(av[i + 1]);
      if ( oldenv )
        {
        env += av[i + 3];
        env += oldenv;
        }
      itksys::SystemTools::PutEnv( env.c_str() );
      i += 4;
      }
    else if ( !skip && strcmp(av[i], "--compare") == 0 )
      {
      if ( i + 2 >= ac )
        {
        usage();
        return 1;
        }
      compareList.push_back( ComparePairType(av[i + 1], av[i + 2]) );
      i += 3;
      }
    else if ( !skip && strcmp(av[i], "--") == 0 )
      {
      skip = true;
      i += 1;
      }
    else if ( !skip && strcmp(av[i], "--help") == 0 )
      {
      usage();
      return 0;
      }
    else
      {
      args.push_back(av[i]);
      i += 1;
      }
    }

  if ( args.empty() )
    {
    usage();
    return 1;
    }

  // a NULL is required at the end of the table
  char **argv = new char *[args.size() + 1];
  for ( i = 0; i < static_cast< int >( args.size() ); i++ )
    {
    argv[i] = args[i];
    }
  argv[args.size()] = NULL;

  itksysProcess *process = itksysProcess_New();
  itksysProcess_SetCommand(process, argv);
  itksysProcess_SetPipeShared(process, itksysProcess_Pipe_STDOUT, true);
  itksysProcess_SetPipeShared(process, itksysProcess_Pipe_STDERR, true);
  itksysProcess_Execute(process);
  itksysProcess_WaitForExit(process, NULL);

  delete[] argv;

  int retCode = itksysProcess_GetExitValue(process);
  if ( retCode != 0 )
    {
    // no need to compare the images: the test has failed
    return retCode;
    }

  // now compare the images
  try
    {
    for ( i = 0; i < static_cast< int >( compareList.size() ); i++ )
      {
      char *testFilename = compareList[i].first;
      char *baselineFilename = compareList[i].second;
      std::cout << "testFilename: " << testFilename << "  baselineFilename: " << baselineFilename << std::endl;

      // Make a list of possible baselines
      std::map< std::string, int >           baselines = RegressionTestBaselines(baselineFilename);
      std::map< std::string, int >::iterator baseline = baselines.begin();
      std::string                            bestBaseline;
      int                                    bestBaselineStatus = itk::NumericTraits< int >::max();
      while ( baseline != baselines.end() )
        {
        baseline->second = RegressionTestImage(testFilename,
                                               ( baseline->first ).c_str(),
                                               0);
        if ( baseline->second < bestBaselineStatus )
          {
          bestBaseline = baseline->first;
          bestBaselineStatus = baseline->second;
          }
        if ( baseline->second == 0 )
          {
          break;
          }
        ++baseline;
        }
      // if the best we can do still has errors, generate the error images
      if ( bestBaselineStatus )
        {
        baseline->second = RegressionTestImage(testFilename,
                                               bestBaseline.c_str(),
                                               1);
        }

      // output the matching baseline
      std::cout << "<DartMeasurement name=\"BaselineImageName\" type=\"text/string\">";
      std::cout << itksys::SystemTools::GetFilenameName(bestBaseline);
      std::cout << "</DartMeasurement>" << std::endl;

      if ( bestBaselineStatus != 0 )
        {
        return bestBaselineStatus;
        }
      }
    }
  catch ( const itk::ExceptionObject & e )
    {
    std::cerr << "ITK test driver caught an ITK exception:\n";
    std::cerr << e.GetFile() << ":" << e.GetLine() << ":\n"
              << e.GetDescription() << "\n";
    return -1;
    }
  catch ( const std::exception & e )
    {
    std::cerr << "ITK test driver caught an exception:\n";
    std::cerr << e.what() << "\n";
    return -1;
    }
  catch ( ... )
    {
    std::cerr << "ITK test driver caught an unknown exception!!!\n";
    return -1;
    }

  return 0;
}
