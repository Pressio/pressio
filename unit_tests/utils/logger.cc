
#include <gtest/gtest.h>
#include "pressio_utils.hpp"

TEST(utils_basic, Logger)
{
  // current choices: terminal, fileAndTerminal, file
  // by default it create mt-safe logger
  pressio::log::initialize(pressio::logto::fileAndTerminal, "log.txt");
  pressio::log::setVerbosity({pressio::log::level::info, pressio::log::level::warn});

  double a = 2.;
  PRESSIOLOG_WARN("this should appear in both console and file");
  PRESSIOLOG_INFO("this message should not appear in the console, only in the file, {:3}", a);
}
