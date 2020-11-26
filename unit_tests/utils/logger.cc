
#include <gtest/gtest.h>
#include "pressio_utils.hpp"

TEST(utils_basic, Logger)
{
  // using namespace pressio;

  // auto a1 = "fr";
  // int b1 = 2;
  // double b2 = 44.4;

  // // auto bg = utils::io::bg_grey();
  // // auto col1 = utils::io::green();
  // // auto reset = utils::io::reset();
  // // utils::io::print_stdout(bg, col1, a1, b1, b2,
  // // 			 reset, "\n");
  // auto logger = pressio::log::create("log.txt");
  // logger.setConsoleLevel(pressio::log::level::warn);
  // logger.setFileLevel(pressio::log::level::info);

  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  console_sink->set_level(spdlog::level::trace);
  //console_sink->set_pattern("[multi_sink_example] [%^%l%$] %v");

  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("multisink.txt", true);
  file_sink->set_level(spdlog::level::info);

  spdlog::sinks_init_list sink_list = { file_sink, console_sink };
  spdlog::logger logger("multi_sink", sink_list.begin(), sink_list.end());
  spdlog::set_default_logger(std::make_shared<spdlog::logger>("multi_sink", spdlog::sinks_init_list({console_sink, file_sink})));

  // spdlog::warn("this should appear in both console and file");
  // spdlog::info("this message should not appear in the console, only in the file");
  SPDLOG_WARN("this should appear in both console and file");
  double a = 2.;
  SPDLOG_INFO("this message should not appear in the console, only in the file {:3}", a);

}
