
#ifndef ODE_TESTING_FUNCTIONAL_SMALL_STEPPERS_CUSTOM_IND_VAR_HPP_
#define ODE_TESTING_FUNCTIONAL_SMALL_STEPPERS_CUSTOM_IND_VAR_HPP_

class MyCustomTime{
  double time_;

public:
  MyCustomTime() = default;
  explicit MyCustomTime(double val) : time_(val){}

  operator double() const{ return time_; }

  MyCustomTime & operator +=(const MyCustomTime & o){
    time_ += o.time_;
    return *this;
  }
};

#endif
