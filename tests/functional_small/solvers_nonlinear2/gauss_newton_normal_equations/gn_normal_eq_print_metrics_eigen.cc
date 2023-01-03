
#include "pressio/solvers.hpp"

struct FakeProblem
{
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;

  state_type createState() const{return state_type(4);}
  residual_type createResidual() const{return residual_type(10);}
  jacobian_type createJacobian() const{return jacobian_type(10,4);}

  void residual(const state_type&, residual_type & res) const
  {
    res.setConstant(1.);
  }

  void jacobian(const state_type &, jacobian_type & jac) const
  {
    jac.setConstant(1.);
  }
};

template<class T>
struct FakeLinS
{
  int count_=0;
  using matrix_type = T;

  template<typename A_t, typename b_t, typename x_t>
  void solve(const A_t &, const b_t &, x_t & x)
  {
    ++count_;
    x.setConstant(1.1);
  }
};

template <typename Out>
void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

int main()
{
  // this is just to test that printing metrics
  // to file works as expected by omitting lables/cols names

  std::string sentinel = "PASSED";
  // Important: use unique log name, because we read the log contents here
  // and don't want any other test to write into the same log file.
  const auto log_file_name = "log_gn_normal_eq_print_metrics_eigen.txt";
  pressio::log::initialize(pressio::logto::fileAndTerminal, log_file_name);
  pressio::log::setVerbosity({pressio::log::level::info, pressio::log::level::debug});

  using namespace pressio;
  using problem_t = FakeProblem;
  using state_t	  = typename problem_t::state_type;
  using mat_t     = typename problem_t::jacobian_type;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(4);
  x(0)=1.0; x(1)=2.; x(2)=3.; x(3)=4.;

  auto solver = pressio::nonlinearsolvers::create_gauss_newton(problem, FakeLinS<hessian_t>{});
  auto criterion = pressio::nonlinearsolvers::Stop::AfterMaxIters;
  solver.setStoppingCriterion(criterion);
  solver.setMaxIterations(3);
  solver.printStrippedMetrics();
  // solve twice
  solver.solve(problem, x);
  solver.solve(problem, x);
  pressio::log::finalize();

  // check file printed
  std::ifstream file(log_file_name);
  if (file.is_open())
  {
    std::string line;
    int lc = 0;
    while (std::getline(file, line))
    {
      std::cout << line << std::endl;
      ++lc;
      if(lc==2){
        // need to start from the end because the beginning of the print 
        // contains the thread ID so that might have different extents 
      	auto s = line.substr(line.size()-28, 28);
        std::cout << s << std::endl;
	      const auto gold = "nonlinsolver: create updater";
	      if (s != gold) sentinel = "FAILED";
      }

      if(lc==3){
      	auto s = line.substr(line.size()-19, 19);
        std::cout << s << std::endl;
      	const auto gold = "nonlinsolver: solve";
      	if (s != gold) sentinel = "FAILED";
      }

      if(lc==4 or lc==5 or lc==6)
      {
      	auto s = line.substr(line.size()-79, 79);
        std::cout << s << std::endl;
      	auto sV = split(s, ' ');
      	if (std::stoi(sV[0]) != (lc-3)) sentinel = "FAILED";
      	if (std::stod(sV[1]) != 3.162278) sentinel = "FAILED";
      	if (std::stod(sV[2]) != 1.) sentinel = "FAILED";
      	if (std::stod(sV[3]) != 20.) sentinel = "FAILED";
      	if (std::stod(sV[4]) != 1.) sentinel = "FAILED";
      	if (std::stod(sV[5]) != 2.2) sentinel = "FAILED";
      	if (std::stod(sV[6]) != 1.) sentinel = "FAILED";
      }

      if(lc==7)
      {
        auto s = line.substr(line.size()-19, 19);
      	const auto gold = "nonlinsolver: solve";
      	if (s != gold) sentinel = "FAILED";
      }
      if(lc==8 or lc==9 or lc==10)
      {
        auto s = line.substr(line.size()-79, 79);
        std::cout << s << std::endl;
      	auto sV = split(s, ' ');
      	if (std::stod(sV[1]) != 3.162278) sentinel = "FAILED";
      	if (std::stod(sV[2]) != 1.) sentinel = "FAILED";
      	if (std::stod(sV[3]) != 20.) sentinel = "FAILED";
      	if (std::stod(sV[4]) != 1.) sentinel = "FAILED";
      	if (std::stod(sV[5]) != 2.2) sentinel = "FAILED";
      	if (std::stod(sV[6]) != 1.) sentinel = "FAILED";
      }

    }
    file.close();
  }

  std::cout << sentinel << std::endl;
  return 0;
}
