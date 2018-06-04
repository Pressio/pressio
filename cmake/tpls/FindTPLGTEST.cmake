
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( GTEST
  #REQUIRED_HEADERS  ...
  REQUIRED_LIBS_NAMES libgmock.dylib libgmock_main.dylib libgtest.dylib libgtest_main.dylib
  MUST_FIND_ALL_HEADERS
  )
