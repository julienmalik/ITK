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
#ifndef __itkGDALImageIO_h
#define __itkGDALImageIO_h

#include "itkImageIOBase.h"
#include <fstream>

#include <gdal.h>

namespace itk
{
//BTX
class GDALDatasetWrapper;
//ETX

/** \class GDALImageIO
 *
 * \brief ImageIO object for reading and writing GDAL images
 *
 * \ingroup IOFilters
 *
 */
class ITK_EXPORT GDALImageIO : public ImageIOBase
{
public:
  /** Standard class typedefs. */
  typedef GDALImageIO          Self;
  typedef ImageIOBase          Superclass;
  typedef SmartPointer< Self > Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GDALImageIO, ImageIOBase);

  /*-------- This part of the interface deals with reading data. ------ */

  /** Determine the file type. Returns true if this ImageIO can read the
   * file specified. */
  virtual bool CanReadFile(const char *);

  /** Set the spacing and diemention information for the set filename. */
  virtual void ReadImageInformation();

  /** Reads the data from disk into the memory buffer provided. */
  virtual void Read(void *buffer);

  /** Reads 3D data from multi-pages GDAL. */
//  virtual void ReadVolume(void *buffer);

  /** Reads 3D data from tiled GDAL. */
//  virtual void ReadTiles(void *buffer);

  /*-------- This part of the interfaces deals with writing data. ----- */

  /** Determine the file type. Returns true if this ImageIO can read the
   * file specified. */
  virtual bool CanWriteFile(const char *);

  /** Determine the file type. Returns true if the ImageIO can stream write the
    specified file */
  virtual bool CanStreamWrite();

  /** Writes the spacing and dimentions of the image.
   * Assumes SetFileName has been called with a valid file name. */
  virtual void WriteImageInformation() {
  }

  /** Writes the data to disk from the memory buffer provided. Make sure
   * that the IORegion has been set properly. */
  virtual void Write(const void *buffer);

protected:
  GDALImageIO();
  ~GDALImageIO();
  void PrintSelf(std::ostream & os, Indent indent) const;

  void InternalWriteImageInformation(const void* buffer);

  void InternalWrite(const void *buffer);

  /** Determine GDAL Driver from the ITK filename */
  std::string FilenameToGdalDriverShortName(std::string name);

  /** Determine the filename to give to GDAL from the ITK filename */
  std::string GetGdalWriteImageFileName(std::string& gdalDriverShortName, std::string filename);

  /** GDAL internal dataset. */
  typedef itk::SmartPointer<GDALDatasetWrapper> GDALDatasetWrapperPointer;
  GDALDatasetWrapperPointer m_Dataset;

  GDALDataType m_PxType;
  int          m_BytePerPixel;
  bool         m_CanStreamWrite;
  bool         m_WriteImageInformationDone;
private:
  GDALImageIO(const Self &);    //purposely not implemented
  void operator=(const Self &); //purposely not implemented

};
} // end namespace itk

#endif // __itkGDALImageIO_h
