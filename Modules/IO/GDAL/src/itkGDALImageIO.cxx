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
#ifdef _MSC_VER
#pragma warning( disable : 4611 )
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <gdal.h>
#include <gdal_priv.h>

#include "itksys/SystemTools.hxx"

#include "itkGDALImageIO.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"
#include "itkTimeProbe.h"
#include "itkMetaDataObject.h"

namespace itk
{

// only two states : the Pointer is Null or GetDataSet() returns a valid dataset
class GDALDatasetWrapper : public itk::LightObject
{
  friend class GDALDriverManagerWrapper;
public:
  typedef GDALDatasetWrapper      Self;
  typedef itk::LightObject        Superclass;
  typedef itk::SmartPointer<Self> Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GDALImageIO, itk::LightObject);

  /** Easy access to the internal GDALDataset object.
   *  Don't close it, it will be automatic */
  GDALDataset* GetDataSet() const
  {
    return m_Dataset;
  }

protected:
  GDALDatasetWrapper()
    : m_Dataset(NULL)
  {
  }

  virtual ~GDALDatasetWrapper()
  {
    if (m_Dataset)
      {
      GDALClose(m_Dataset);
      }
  }

private:
  GDALDataset* m_Dataset;
}; // end of GDALDatasetWrapper

// Wraps the GdalDriverManager so that GDALAllRegister is called automatically
class GDALDriverManagerWrapper
{
public:
  // GetInstance returns a reference to a GDALDriverManagerWrapper
  // This is the only entry point to interact with this class
  static GDALDriverManagerWrapper& GetInstance()
  {

    // Declare a static method variable of type GDALDriverManagerWrapper
    // so that it is constructed and initialized only on the first call
    // to GetInstance(), and so try to avoid static initialization order
    // problems

    static GDALDriverManagerWrapper theUniqueInstance;

    return theUniqueInstance;
  }

  // Open the file for reading and returns a smart dataset pointer
  GDALDatasetWrapper::Pointer Open( std::string filename ) const
  {
    GDALDatasetWrapper::Pointer datasetWrapper;
    GDALDatasetH                dataset = GDALOpen(filename.c_str(), GA_ReadOnly);

    if (dataset != NULL)
      {
      datasetWrapper = GDALDatasetWrapper::New();
      datasetWrapper->m_Dataset = static_cast<GDALDataset*>(dataset);
      }
    return datasetWrapper;
  }

  // Open the new  file for writing and returns a smart dataset pointer
  GDALDatasetWrapper::Pointer Create( std::string driverShortName, std::string filename,
                                      int nXSize, int nYSize, int nBands,
                                      GDALDataType eType, char ** papszOptions ) const
  {
    GDALDatasetWrapper::Pointer datasetWrapper;

    GDALDriver* driver = GetDriverByName( driverShortName );

    if(driver != NULL)
      {
      GDALDataset* dataset = driver->Create(filename.c_str(),
                                            nXSize, nYSize,
                                            nBands, eType,
                                            papszOptions );

      if (dataset != NULL)
        {
        datasetWrapper = GDALDatasetWrapper::New();
        datasetWrapper->m_Dataset = dataset;
        }
      }
    return datasetWrapper;
  }

  GDALDriver* GetDriverByName( std::string driverShortName ) const
  {
    return GetGDALDriverManager()->GetDriverByName(driverShortName.c_str() );
  }

private:
  // private constructor so that this class is allocated only inside GetInstance
  GDALDriverManagerWrapper()
  {
    GDALAllRegister();
  }

  ~GDALDriverManagerWrapper()
  {
    GDALDestroyDriverManager();
  }

}; // end of GDALDriverManagerWrapper

GDALImageIO::GDALImageIO()
{
  // By default set number of dimensions to two.
  this->SetNumberOfDimensions(2);

  // By default set pixel type to scalar.
  m_PixelType = SCALAR;

  // By default set component type to unsigned char
  m_ComponentType = UCHAR;

  // Set default spacing to one
  m_Spacing[0] = 1.0;
  m_Spacing[1] = 1.0;
  // Set default origin to zero
  m_Origin[0] = 0.0;
  m_Origin[1] = 0.0;

  m_CanStreamWrite = false;
  m_WriteImageInformationDone = false;
}

GDALImageIO::~GDALImageIO()
{
}

std::string GDALImageIO::FilenameToGdalDriverShortName(std::string name)
{
  std::string extension;
  std::string gdalDriverShortName;

  // Get extension in lowercase
  extension = itksys::SystemTools::LowerCase( itksys::SystemTools::GetFilenameLastExtension(name) );

  if      ( extension == ".tif" || extension == ".tiff" )
    gdalDriverShortName = "GTiff";
  else if ( extension == ".hdr" )
    gdalDriverShortName = "ENVI";
  else if ( extension == ".img" )
    gdalDriverShortName = "HFA";
  else if ( extension == ".ntf" )
    gdalDriverShortName = "NITF";
  else if ( extension == ".png" )
    gdalDriverShortName="PNG";
  else if ( extension == ".bmp" )
    gdalDriverShortName="BMP";
  else if ( extension == ".j2k" )
    gdalDriverShortName="JPEG2000";
  else if ( extension == ".jpg" || extension==".jpeg" )
    gdalDriverShortName="JPEG";
  else if ( extension == ".pix" )
    gdalDriverShortName="PCIDSK";
  else if ( extension == ".ras" )
    gdalDriverShortName="GTiff";
  else
    gdalDriverShortName = "NOT-FOUND";

  return gdalDriverShortName;
}

std::string GDALImageIO::GetGdalWriteImageFileName(std::string& gdalDriverShortName, std::string filename)
{
  std::string gdalFileName;

  gdalFileName = filename;
  // Suppress hdr extension for ENVI format
  if (gdalDriverShortName == "ENVI")
    {
    gdalFileName = itksys::SystemTools::GetFilenameWithoutLastExtension(filename);
    }
  return gdalFileName;
}

bool GDALImageIO::CanReadFile(const char *file)
{
  // First check the extension
  if (file == NULL)
    {
    itkDebugMacro(<< "No filename specified.");
    return false;
    }

  // Try to open the file to see if GDAL supports it
  m_Dataset = GDALDriverManagerWrapper::GetInstance().Open(file);

  // If we get a valid pointer, it means we can read the file
  return m_Dataset.IsNotNull();
}

void GDALImageIO::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

// Read image with GDAL
void GDALImageIO::Read(void* buffer)
{
  // Convert buffer from void * to unsigned char *
  unsigned char *p = static_cast<unsigned char *>(buffer);

  // Check if conversion succeed
  if (p == NULL)
    {
    itkExceptionMacro(<< "GDAL : Null pointer given for reading");
    return;
    }

  // Get nb. of lines and columns of the region to read
  int nbLines     = this->GetIORegion().GetSize()[1];
  int nbColumns   = this->GetIORegion().GetSize()[0];
  int firstLine   = this->GetIORegion().GetIndex()[1]; // [1... ]
  int firstColumn = this->GetIORegion().GetIndex()[0]; // [1... ]

  GDALDataset* dataset = m_Dataset->GetDataSet();

  int pixelOffset = m_BytePerPixel * this->GetNumberOfComponents();
  int lineOffset  = m_BytePerPixel * this->GetNumberOfComponents() * nbColumns;
  int bandOffset  = m_BytePerPixel;
  int nbBands     = this->GetNumberOfComponents();

  itkDebugMacro(<< "GDAL RasterIO : \n"
                << " indX = " << firstColumn << "\n"
                << " indY = " << firstLine << "\n"
                << " sizeX = " << nbColumns << "\n"
                << " sizeY = " << nbLines << "\n"
                << " GDAL Data Type = " << GDALGetDataTypeName(m_PxType) << "\n"
                << " pixelOffset = " << pixelOffset << "\n"
                << " lineOffset = " << lineOffset << "\n"
                << " bandOffset = " << bandOffset)

  itk::TimeProbe chrono;
  chrono.Start();
  CPLErr status = dataset->RasterIO(GF_Read,
                                    firstColumn,
                                    firstLine,
                                    nbColumns,
                                    nbLines,
                                    p,
                                    nbColumns,
                                    nbLines,
                                    m_PxType,
                                    nbBands,
                                    NULL,  // We want to read all bands
                                    pixelOffset,
                                    lineOffset,
                                    bandOffset);
  chrono.Stop();
  itkDebugMacro(<< "RasterIO Read took " << chrono.GetTotal() << " sec")

  // Check if gdal call succeed
  if (status == CE_Failure)
    {
    itkExceptionMacro(<< "Error while reading image (GDAL format) " << m_FileName.c_str() << ".")
    }
}

void GDALImageIO::ReadImageInformation()
{
  if ( !m_Dataset )
    m_Dataset = GDALDriverManagerWrapper::GetInstance().Open(this->GetFileName() );

  GDALDataset* dataset = m_Dataset->GetDataSet();

  // Get image dimensions
  if ( dataset->GetRasterXSize() == 0 || dataset->GetRasterYSize() == 0 )
    {
    itkExceptionMacro(<< "Dimension is undefined.");
    }

  // Set image dimensions into IO
  m_Dimensions[0] = dataset->GetRasterXSize();
  m_Dimensions[1] = dataset->GetRasterYSize();

  this->SetNumberOfComponents(dataset->GetRasterCount() );

  itkDebugMacro(<< "Input file dimension: " << m_Dimensions[0] << ", " << m_Dimensions[1]);
  itkDebugMacro(<< "Number of bands inside input file: " << this->GetNumberOfComponents() );

  // Set the number of dimensions (verify for the dim )
  this->SetNumberOfDimensions(2);

  // Automatically set the Type to Binary for GDAL data
  this->SetFileTypeToBinary();

  // Get Data Type
  // Consider only the data type given by the first band
  // Maybe should changed (to check)
  m_PxType = dataset->GetRasterBand(1)->GetRasterDataType();
  itkDebugMacro(<< "PixelType inside input file: "<< GDALGetDataTypeName(m_PxType) );
  if (m_PxType == GDT_Byte)
    {
    SetComponentType(UCHAR);
    }
  else if (m_PxType == GDT_UInt16)
    {
    SetComponentType(USHORT);
    }
  else if (m_PxType == GDT_Int16)
    {
    SetComponentType(SHORT);
    }
  else if (m_PxType == GDT_UInt32)
    {
    SetComponentType(UINT);
    }
  else if (m_PxType == GDT_Int32)
    {
    SetComponentType(INT);
    }
  else if (m_PxType == GDT_Float32)
    {
    SetComponentType(FLOAT);
    }
  else if (m_PxType == GDT_Float64)
    {
    SetComponentType(DOUBLE);
    }
  /*
  else if (m_PxType == GDT_CInt16)
    {
    SetComponentType(CSHORT);
    }
  else if (m_PxType == GDT_CInt32)
    {
    SetComponentType(CINT);
    }
  else if (m_PxType == GDT_CFloat32)
    {
    SetComponentType(CFLOAT);
    }
  else if (m_PxType == GDT_CFloat64)
    {
    SetComponentType(CDOUBLE);
    }
    */
  else
    {
    itkExceptionMacro(<< "Pixel type unknown");
    }

  if (this->GetComponentType() == CHAR)
    {
    m_BytePerPixel = 1;
    }
  else if (this->GetComponentType() == UCHAR)
    {
    m_BytePerPixel = 1;
    }
  else if (this->GetComponentType() == USHORT)
    {
    m_BytePerPixel = 2;
    }
  else if (this->GetComponentType() == SHORT)
    {
    m_BytePerPixel = 2;
    }
  else if (this->GetComponentType() == INT)
    {
    m_BytePerPixel = 4;
    }
  else if (this->GetComponentType() == UINT)
    {
    m_BytePerPixel = 4;
    }
  else if (this->GetComponentType() == FLOAT)
    {
    m_BytePerPixel = 4;
    }
  else if (this->GetComponentType() == DOUBLE)
    {
    m_BytePerPixel = 8;
    }
  /*
  else if (this->GetComponentType() == CSHORT)
    {
    m_BytePerPixel = sizeof(std::complex<short>);
    }
  else if (this->GetComponentType() == CINT)
    {
    m_BytePerPixel = sizeof(std::complex<int>);
    }
  else if (this->GetComponentType() == CFLOAT)
    {
    m_BytePerPixel = sizeof(std::complex<float>);
    }
  else if (this->GetComponentType() == CDOUBLE)
    {
      m_BytePerPixel = sizeof(std::complex<double>);
    }
  */
  else
    {
    itkExceptionMacro(<< "Component type unknown");
    }

  if (this->GetNumberOfComponents() == 1)
    {
    this->SetPixelType(SCALAR);
    }
  else
    {
    this->SetPixelType(VECTOR);
    }

  /*----------------------------------------------------------------------*/
  /*-------------------------- METADATA ----------------------------------*/
  /*----------------------------------------------------------------------*/
  itk::MetaDataDictionary& dict = this->GetMetaDataDictionary();

  /* -------------------------------------------------------------------- */
  /*  Get Spacing  and origin                                             */
  /* -------------------------------------------------------------------- */
  // Default Spacing
  m_Spacing[0] = 1;
  m_Spacing[1] = 1;
  if (m_NumberOfDimensions == 3) m_Spacing[2] = 1;
  double geoTransform[6];
  if (dataset->GetGeoTransform(geoTransform) == CE_None)
    {
    m_Origin[0]  = geoTransform[0];
    m_Origin[1]  = geoTransform[3];
    m_Spacing[0] = geoTransform[1];
    m_Spacing[1] = geoTransform[5];
    }

  /* -------------------------------------------------------------------- */
  /*      Report driver info.                                             */
  /* -------------------------------------------------------------------- */
  GDALDriverH hDriver;

  hDriver = dataset->GetDriver();

  std::string driverShortName =  static_cast<std::string>(GDALGetDriverShortName(hDriver) );
  std::string driverLongName  =  static_cast<std::string>(GDALGetDriverLongName(hDriver) );

  itk::EncapsulateMetaData<std::string>(dict, "DriverShortName", driverShortName);
  itk::EncapsulateMetaData<std::string>(dict, "DriverLongName",  driverLongName);

  /* -------------------------------------------------------------------- */
  /*      Report metadata.                                                */
  /* -------------------------------------------------------------------- */

  char** papszMetadata;
  papszMetadata = dataset->GetMetadata(NULL);
  if (CSLCount(papszMetadata) > 0)
    {
    std::string key;

    for (int cpt = 0; papszMetadata[cpt] != NULL; cpt++)
      {
      ::std::ostringstream ost;
      ost << "Metadata_" << cpt;
      key = ost.str();

      itk::EncapsulateMetaData<std::string>(dict, key,
                                            static_cast<std::string>(papszMetadata[cpt]) );
      }
    }

}

bool GDALImageIO::CanWriteFile(const char* name)
{
  m_FileName = name;

  // First check the filename
  if (name == NULL)
    {
    itkDebugMacro(<< "No filename specified.");
    return false;
    }

  // Get the GDAL format ID from the name
  std::string gdalDriverShortName = FilenameToGdalDriverShortName(name);
  if (gdalDriverShortName == "NOT-FOUND")
    {
    itkDebugMacro( "No suitable GDAL driver found for " << name )
    return false;
    }
  itkDebugMacro( "Found GDAL driver : " << gdalDriverShortName << std::endl )

  // Check the driver for support of Create or at least CreateCopy
  GDALDriver* driver = GDALDriverManagerWrapper::GetInstance().GetDriverByName(gdalDriverShortName);
  if ( GDALGetMetadataItem( driver, GDAL_DCAP_CREATE, NULL ) == NULL
       && GDALGetMetadataItem( driver, GDAL_DCAP_CREATECOPY, NULL ) == NULL )
    {
    itkDebugMacro(<< "The driver " << GDALGetDriverShortName(driver) << " does not support writing");
    return false;
    }
  return true;
}

bool GDALImageIO::CanStreamWrite()
{
  // Get the GDAL format ID from the name
  std::string gdalDriverShortName = FilenameToGdalDriverShortName(m_FileName);
  GDALDriver* driver = GDALDriverManagerWrapper::GetInstance().GetDriverByName(gdalDriverShortName);

  if (driver == NULL)
    {
    itkDebugMacro(<< "Unable to instantiate driver " << gdalDriverShortName);
    m_CanStreamWrite = false;
    }

  if ( GDALGetMetadataItem( driver, GDAL_DCAP_CREATE, NULL ) != NULL )
    {
    m_CanStreamWrite = true;
    }
  else
    {
    m_CanStreamWrite = false;
    }
  return m_CanStreamWrite;
}

void GDALImageIO::InternalWriteImageInformation(const void* buffer)
{
  // Check if we have to write the image information
  if (m_WriteImageInformationDone == true)
    {
    return;
    }

  char **     papszOptions = NULL;
  std::string driverShortName;

  if ( m_Dimensions[0] == 0 || m_Dimensions[1] == 0 )
    {
    itkExceptionMacro(<< "Null Dimension");
    }

  if (this->GetComponentType() == CHAR)
    {
    m_BytePerPixel = 1;
    m_PxType = GDT_Byte;
    }
  else if (this->GetComponentType() == UCHAR)
    {
    m_BytePerPixel = 1;
    m_PxType = GDT_Byte;
    }
  else if (this->GetComponentType() == USHORT)
    {
    m_BytePerPixel = 2;
    m_PxType = GDT_UInt16;
    }
  else if (this->GetComponentType() == SHORT)
    {
    m_BytePerPixel = 2;
    m_PxType = GDT_Int16;
    }
  else if (this->GetComponentType() == INT)
    {
    m_BytePerPixel = 4;
    m_PxType = GDT_Int32;
    }
  else if (this->GetComponentType() == UINT)
    {
    m_BytePerPixel = 4;
    m_PxType = GDT_UInt32;
    }
  else if (this->GetComponentType() == FLOAT)
    {
    m_BytePerPixel = 4;
    m_PxType = GDT_Float32;
    }
  else if (this->GetComponentType() == DOUBLE)
    {
    m_BytePerPixel = 8;
    m_PxType = GDT_Float64;
    }
  else
    {
    m_BytePerPixel = 1;
    m_PxType = GDT_Byte;
    }

  // Automatically set the Type to Binary for GDAL data
  this->SetFileTypeToBinary();

  driverShortName = FilenameToGdalDriverShortName(m_FileName);
  if (driverShortName == "NOT-FOUND")
    {
    itkExceptionMacro(
      << "GDAL Writing failed : the image file name '" << m_FileName.c_str() << "' is not recognized by GDAL.");
    }

  if (m_CanStreamWrite)
    {
    // Force tile mode for TIFF format. Tile mode is a lot more
    // efficient when writing huge tiffs
    if( driverShortName.compare("GTiff") == 0 )
      {
      itkDebugMacro(<< "Enabling TIFF Tiled mode")
      papszOptions = CSLAddNameValue( papszOptions, "TILED", "YES" );

      // Use a fixed tile size
      // Take as reference is a 256*256 short int 4 bands tile
      const unsigned int ReferenceTileSizeInBytes = 256 * 256 * 4 * 2;

      unsigned int nbPixelPerTile = ReferenceTileSizeInBytes / m_BytePerPixel / this->GetNumberOfComponents();
      unsigned int tileDimension = static_cast<unsigned int>( vcl_sqrt(static_cast<float>(nbPixelPerTile) ) );

      // align the tile dimension to the next multiple of 16 (needed by TIFF
      // spec)
      tileDimension = ( tileDimension + 15 ) / 16 * 16;

      itkDebugMacro(<< "Tile dimension : " << tileDimension << " * " << tileDimension)

      std::ostringstream oss;
      oss << tileDimension;
      papszOptions = CSLAddNameValue( papszOptions, "BLOCKXSIZE", oss.str().c_str() );
      papszOptions = CSLAddNameValue( papszOptions, "BLOCKYSIZE", oss.str().c_str() );
      }

    m_Dataset = GDALDriverManagerWrapper::GetInstance().Create(
        driverShortName,
        GetGdalWriteImageFileName(driverShortName, m_FileName),
        m_Dimensions[0], m_Dimensions[1],
        this->GetNumberOfComponents(), m_PxType,
        papszOptions);
    }
  else
    {
    // buffer casted in unsigned long cause under Win32 the adress
    // don't begin with 0x, the adress in not interpreted as
    // hexadecimal but alpha numeric value, then the conversion to
    // integer make us pointing to an non allowed memory block => Crash.
    std::ostringstream stream;
    stream << "MEM:::"
           <<  "DATAPOINTER=" << (unsigned long)(buffer) << ","
           <<  "PIXELS=" << m_Dimensions[0] << ","
           <<  "LINES=" << m_Dimensions[1] << ","
           <<  "BANDS=" << this->GetNumberOfComponents() << ","
           <<  "DATATYPE=" << GDALGetDataTypeName(m_PxType) << ","
           <<  "PIXELOFFSET=" << m_BytePerPixel * this->GetNumberOfComponents() << ","
           <<  "LINEOFFSET=" << m_BytePerPixel * this->GetNumberOfComponents() * m_Dimensions[0] << ","
           <<  "BANDOFFSET=" << m_BytePerPixel;

    m_Dataset = GDALDriverManagerWrapper::GetInstance().Open(stream.str() );
    }

  if (m_Dataset.IsNull() )
    {
    itkExceptionMacro(
      << "GDAL Writing failed : Impossible to create the image file name '" << m_FileName << "'.");
    }

  /*----------------------------------------------------------------------*/
  /*-------------------------- METADATA ----------------------------------*/
  /*----------------------------------------------------------------------*/
  GDALDataset* dataset = m_Dataset->GetDataSet();

  /* -------------------------------------------------------------------- */
  /*  Set the six coefficients of affine geotransform                     */
  /* -------------------------------------------------------------------- */
  itk::VariableLengthVector<double> geoTransform(6);
  /// Reporting origin and spacing
  geoTransform[0] = m_Origin[0];
  geoTransform[3] = m_Origin[1];
  geoTransform[1] = m_Spacing[0];
  geoTransform[5] = m_Spacing[1];
  // FIXME: Here component 1 and 4 should be replaced by the orientation
  // parameters
  geoTransform[2] = 0.;
  geoTransform[4] = 0.;
  dataset->SetGeoTransform(const_cast<double*>(geoTransform.GetDataPointer() ) );

  m_WriteImageInformationDone = true;
}

void GDALImageIO::Write(const void* buffer)
{
  this->InternalWriteImageInformation(buffer);

  // Compute offset and size
  unsigned int nbLines = this->GetIORegion().GetSize()[1];
  unsigned int nbColumns = this->GetIORegion().GetSize()[0];
  int          firstLine = this->GetIORegion().GetIndex()[1];   // [1... ]
  int          firstColumn = this->GetIORegion().GetIndex()[0]; // [1... ]

  // If driver supports streaming
  if (m_CanStreamWrite)
    {

    itkDebugMacro(<< "RasterIO Write requested region : " << this->GetIORegion() <<
                  "\n, lNbColumns =" << nbColumns <<
                  "\n, lNbLines =" << nbLines <<
                  "\n, m_PxType =" << GDALGetDataTypeName(m_PxType) <<
                  "\n, m_NbBands =" << this->GetNumberOfComponents() <<
                  "\n, m_BytePerPixel ="<< m_BytePerPixel <<
                  "\n, Pixel offset =" << m_BytePerPixel * this->GetNumberOfComponents() <<
                  "\n, Line offset =" << m_BytePerPixel * this->GetNumberOfComponents() * nbColumns <<
                  "\n, Band offset =" <<  m_BytePerPixel)
    itk::TimeProbe chrono;
    chrono.Start();
    CPLErr status = m_Dataset->GetDataSet()->RasterIO(GF_Write,
                                                      firstColumn,
                                                      firstLine,
                                                      nbColumns,
                                                      nbLines,
                                                      const_cast<void *>(buffer),
                                                      nbColumns,
                                                      nbLines,
                                                      m_PxType,
                                                      this->GetNumberOfComponents(),
                                                      NULL,
                                                      m_BytePerPixel * this->GetNumberOfComponents(),
                                                      m_BytePerPixel * this->GetNumberOfComponents() * nbColumns,
                                                      m_BytePerPixel);
    chrono.Stop();
    itkDebugMacro(<< "RasterIO Write took " << chrono.GetTotal() << " sec")

    // Check if writing succeed
    if (status == CE_Failure)
      {
      itkExceptionMacro(<< "Error while writing image (GDAL format) " << m_FileName.c_str() << ".");
      }
    // Flush dataset cache
    m_Dataset->GetDataSet()->FlushCache();
    }
  else
    {
    // We only wrote data to the memory dataset
    // Now write it to the real file with CreateCopy()
    std::string gdalDriverShortName = FilenameToGdalDriverShortName(m_FileName);
    std::string realFileName = GetGdalWriteImageFileName(gdalDriverShortName, m_FileName);

    GDALDriver* driver = GDALDriverManagerWrapper::GetInstance().GetDriverByName(gdalDriverShortName);
    if (driver == NULL)
      {
      itkExceptionMacro(<< "Unable to instantiate driver " << gdalDriverShortName << " to write " << m_FileName);
      }

    // If JPEG, set the JPEG compression quality to 95.
    char * option[2];
    option[0] = NULL;
    option[1] = NULL;
    // If JPEG, set the image quality
    if( gdalDriverShortName.compare("JPEG") == 0 )
      {
      option[0] = const_cast<char *>("QUALITY=95");
      }

    GDALDataset* hOutputDS = driver->CreateCopy( realFileName.c_str(), m_Dataset->GetDataSet(), FALSE,
                                                 option, NULL, NULL );
    GDALClose(hOutputDS);
    }
}

} // end namespace itk
