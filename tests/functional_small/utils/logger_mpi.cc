
#include <gtest/gtest.h>
#include "pressio/utils.hpp"

TEST(utils_basic, loggerMpi)
{
  // current choices: terminal, fileAndTerminal, file
  // by default it create mt-safe logger
  pressio::log::initialize(pressio::logto::fileAndTerminal, "log.txt");
  pressio::log::setVerbosity({pressio::log::level::info, pressio::log::level::warn});

  double a = 2.;
  PRESSIOLOG_WARN("this message should appear in both console and file");
  PRESSIOLOG_INFO("this message should not appear in the console, only in the file {:3}", a);
  pressio::log::finalize();
}
