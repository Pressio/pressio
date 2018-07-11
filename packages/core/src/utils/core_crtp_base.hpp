
#ifndef CORE_CRPT_BASE_HPP_
#define CORE_CRTP_BASE_HPP_ 

namespace core{
  
template <typename T,
	  template<typename, typename...> class crtpType,
	  typename ... Args>
struct crtpBase{
  T & underlying() {
    return static_cast<T&>(*this);
  }
  T const & underlying() const {
    return static_cast<T const&>(*this);
  }
private:
  crtpBase(){}
  friend crtpType<T, Args...>;

};//end class
}//end namespace
#endif
