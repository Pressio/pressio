#include <pressio/utils.hpp>

int main() {
  pressio::log::initialize(pressio::logto::terminal, "log.txt");
  pressio::log::setVerbosity({pressio::log::level::info});

  double a = 1.;
  PRESSIOLOG_INFO("pressio log message, {:3}", a);

  pressio::log::finalize();
}
