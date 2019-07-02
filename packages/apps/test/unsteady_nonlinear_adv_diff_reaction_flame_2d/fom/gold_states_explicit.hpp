
#ifndef PRESSIO_APPS_TEST_UNSTEADY_NONLIN_ADVDIFFREACTION_FLAME_2D_GOLD_EXPLICIT_HPP_
#define PRESSIO_APPS_TEST_UNSTEADY_NONLIN_ADVDIFFREACTION_FLAME_2D_GOLD_EXPLICIT_HPP_

#include "ODE_ALL"

namespace pressio { namespace apps{ namespace test{

template <ode::ExplicitEnum>
struct NonLinAdvDiffReacFlame2dExpGoldStates;

template <>
struct NonLinAdvDiffReacFlame2dExpGoldStates<ode::ExplicitEnum::RungeKutta4>{
  using result_t = std::vector<double>;

  static result_t get(int Nx, int Ny, double dt, double final_t){
    if(Nx==12 and Ny==6 and
       std::abs(dt-0.00001)<1e-6 and
       std::abs(final_t-0.011)<1e-6)
      {

      return{
          329.941462339757265,
            0.000370644911049,
            0.002988251801415,
            0.001728233605795,
          430.478026085137344,
            0.001383342390631,
            0.011169686894078,
            0.007962847352763,
          464.031694423233034,
            0.001929335368813,
            0.015562204269555,
            0.009656391891988,
          428.620343000864978,
            0.001770764537869,
            0.014263531633000,
            0.007091581862210,
          375.387246909466285,
            0.001211429938769,
            0.009746803249789,
            0.003833460025911,
          335.812028756429982,
            0.000658639991993,
            0.005294589277580,
            0.001666236728536,
          314.408007516615783,
            0.000296984480004,
            0.002385803111803,
            0.000610797055694,
          305.045747915681659,
            0.000114456501432,
            0.000919024188081,
            0.000194447718213,
          301.567481104653211,
            0.000038543032441,
            0.000309362005157,
            0.000054845904641,
          300.437989139856541,
            0.000011530909718,
            0.000092524055598,
            0.000013908243864,
          300.111798958555482,
            0.000003120999330,
            0.000025036868782,
            0.000003219386719,
          300.018354097049212,
            0.000000539089650,
            0.000004323757463,
            0.000000478798443,
          443.229699421897578,
            0.001753294533912,
            0.014137020164736,
            0.008304032754767,
          657.444397487657056,
            0.003208164449586,
            0.025953149525909,
            0.022896738220074,
          669.425811003251738,
            0.003514756245978,
            0.028413543618345,
            0.023293665684517,
          558.322620626983507,
            0.002953158377999,
            0.023826943183992,
            0.015365863661307,
          440.183897608526195,
            0.001947029229746,
            0.015682236787690,
            0.007697379176383,
          363.132151224094287,
            0.001037147556046,
            0.008343317900125,
            0.003168119925219,
          324.452547427052139,
            0.000461274826759,
            0.003707473736923,
            0.001116202314268,
          308.327365190685043,
            0.000175978579184,
            0.001413523013888,
            0.000344958248103,
          302.532684991159385,
            0.000058794920450,
            0.000472038582357,
            0.000095099720638,
          300.696117722027566,
            0.000017478605105,
            0.000140277208622,
            0.000023683773750,
          300.175358006982663,
            0.000004706054007,
            0.000037758352616,
            0.000005401995338,
          300.028479282158230,
            0.000000809265016,
            0.000006491522719,
            0.000000793600346,
         1643.965455904206237,
            0.016474080451866,
            0.132830810857490,
            0.077877401657651,
         1717.440348107389809,
            0.009707395333059,
            0.078830676076799,
            0.096408554388075,
         1374.630235062671773,
            0.007186103690104,
            0.058378744157931,
            0.073415068250086,
          933.432429008688814,
            0.005393592418660,
            0.043661727285669,
            0.041118568638780,
          606.951698849492686,
            0.003448842086172,
            0.027830980071721,
            0.018370633954474,
          428.164795922180929,
            0.001807574651346,
            0.014557237169143,
            0.006986263379627,
          347.159283388336860,
            0.000793752729406,
            0.006384296296677,
            0.002331173437499,
          315.488753831291945,
            0.000299692030773,
            0.002408403642437,
            0.000693043412876,
          304.586881250013732,
            0.000099281891784,
            0.000797368246652,
            0.000185636446302,
          301.235433217289938,
            0.000029309176367,
            0.000235286241172,
            0.000045217814438,
          300.306295875830529,
            0.000007845286750,
            0.000062958034164,
            0.000010133197705,
          300.049108569175758,
            0.000001342370734,
            0.000010769504592,
            0.000001467294861,
         1643.965455904206237,
            0.016474080451866,
            0.132830810857490,
            0.077877401657651,
         1717.440348107389809,
            0.009707395333059,
            0.078830676076799,
            0.096408554388075,
         1374.630235062671773,
            0.007186103690104,
            0.058378744157931,
            0.073415068250086,
          933.432429008688814,
            0.005393592418660,
            0.043661727285669,
            0.041118568638780,
          606.951698849492686,
            0.003448842086172,
            0.027830980071721,
            0.018370633954474,
          428.164795922180929,
            0.001807574651346,
            0.014557237169143,
            0.006986263379627,
          347.159283388336860,
            0.000793752729406,
            0.006384296296677,
            0.002331173437499,
          315.488753831291945,
            0.000299692030773,
            0.002408403642437,
            0.000693043412876,
          304.586881250013732,
            0.000099281891784,
            0.000797368246652,
            0.000185636446302,
          301.235433217289938,
            0.000029309176367,
            0.000235286241172,
            0.000045217814438,
          300.306295875830529,
            0.000007845286750,
            0.000062958034164,
            0.000010133197705,
          300.049108569175758,
            0.000001342370734,
            0.000010769504592,
            0.000001467294861,
          443.229699421897578,
            0.001753294533912,
            0.014137020164736,
            0.008304032754767,
          657.444397487657056,
            0.003208164449586,
            0.025953149525909,
            0.022896738220074,
          669.425811003251738,
            0.003514756245978,
            0.028413543618345,
            0.023293665684517,
          558.322620626983507,
            0.002953158377999,
            0.023826943183992,
            0.015365863661307,
          440.183897608526195,
            0.001947029229746,
            0.015682236787690,
            0.007697379176383,
          363.132151224094287,
            0.001037147556046,
            0.008343317900125,
            0.003168119925219,
          324.452547427052139,
            0.000461274826759,
            0.003707473736923,
            0.001116202314268,
          308.327365190685043,
            0.000175978579184,
            0.001413523013888,
            0.000344958248103,
          302.532684991159385,
            0.000058794920450,
            0.000472038582357,
            0.000095099720638,
          300.696117722027566,
            0.000017478605105,
            0.000140277208622,
            0.000023683773750,
          300.175358006982663,
            0.000004706054007,
            0.000037758352616,
            0.000005401995338,
          300.028479282158230,
            0.000000809265016,
            0.000006491522719,
            0.000000793600346,
          329.941462339757265,
            0.000370644911049,
            0.002988251801415,
            0.001728233605795,
          430.478026085137344,
            0.001383342390631,
            0.011169686894078,
            0.007962847352763,
          464.031694423233034,
            0.001929335368813,
            0.015562204269555,
            0.009656391891988,
          428.620343000864978,
            0.001770764537869,
            0.014263531633000,
            0.007091581862210,
          375.387246909466285,
            0.001211429938769,
            0.009746803249789,
            0.003833460025911,
          335.812028756429982,
            0.000658639991993,
            0.005294589277580,
            0.001666236728536,
          314.408007516615783,
            0.000296984480004,
            0.002385803111803,
            0.000610797055694,
          305.045747915681659,
            0.000114456501432,
            0.000919024188081,
            0.000194447718213,
          301.567481104653211,
            0.000038543032441,
            0.000309362005157,
            0.000054845904641,
          300.437989139856541,
            0.000011530909718,
            0.000092524055598,
            0.000013908243864,
          300.111798958555482,
            0.000003120999330,
            0.000025036868782,
            0.000003219386719,
          300.018354097049212,
            0.000000539089650,
            0.000004323757463,
            0.000000478798443
      };//end return
    }else{
      std::cout << "returning empty true solution" << std::endl;
      return {};
    }
  }//get
};//struct

}}}//end namespace pressio::apps::test
#endif
