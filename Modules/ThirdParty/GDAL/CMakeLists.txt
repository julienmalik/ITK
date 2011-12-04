project(ITKGDAL)
set(ITKGDAL_THIRD_PARTY 1)

option(ITK_USE_SYSTEM_GDAL "Use system-installed gdal" ON)
mark_as_advanced(ITK_USE_SYSTEM_GDAL)

if(ITK_USE_SYSTEM_GDAL)
  find_package(GDAL REQUIRED)
#  set(ITK-GDAL_INCLUDE_DIRS
#    ${ITK-GDAL_BINARY_DIR}/src
#    )
  set(ITKGDAL_SYSTEM_INCLUDE_DIRS
    ${GDAL_INCLUDE_DIRS}
    )
  set(ITKGDAL_LIBRARIES "${GDAL_LIBRARIES}")
  set(ITKGDAL_NO_SRC 1)
else()
  message(FATAL_ERROR "ITK_USE_SYSTEM_GDAL=OFF is not supported. Please install gdal or deactivate the ITKGDAL module")
endif()

itk_module_impl()
