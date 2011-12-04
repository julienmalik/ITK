set(ITKIOGDAL_SRC
itkGDALImageIO.cxx
itkGDALImageIOFactory.cxx
)

add_library(ITKIOGDAL ${ITKIOGDAL_SRC})
target_link_libraries(ITKIOGDAL ${ITKGDAL_LIBRARIES} ${ITKIOImageBase_LIBRARIES})
itk_module_target(ITKIOGDAL)
