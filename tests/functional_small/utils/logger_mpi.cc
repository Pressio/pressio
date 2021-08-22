
#include <mpi.h>
#include <gtest/gtest.h>
#include "pressio/utils.hpp"

TEST(utils_basic, loggerMpi)
{
  int rank = {};
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // current choices: terminal, fileAndTerminal, file
  // by default it create mt-safe logger

  // const auto log_file = "my_log" + std::to_string(rank) + ".txt";
  pressio::log::initialize(pressio::logto::fileAndTerminal, "log.txt");
  pressio::log::setVerbosity({pressio::log::level::info, pressio::log::level::warn});

  double a = 2.;
  PRESSIOLOG_WARN("message from rank = {:2} in both console and file", rank);
  PRESSIOLOG_INFO("this message should not appear in the console, only in the file {:3}", a);

  pressio::log::finalize();
}
