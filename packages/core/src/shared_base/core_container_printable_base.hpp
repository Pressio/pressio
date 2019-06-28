
#ifndef CORE_SHARED_BASE_CONTAINER_PRINTABLE_BASE_HPP_
#define CORE_SHARED_BASE_CONTAINER_PRINTABLE_BASE_HPP_

#include "../core_ConfigDefs.hpp"

namespace rompp{ namespace core{

/*
 * number of elements to prints defaults = -1 , indicating to print all
 *
 * the char c = 'd', 'f'
 * 'd' = default, prints as they are, so for vectors prints vertically
 * 'f' = flatten, for vectors prints horizontally,
 *		  for matrices it falttens it
 */

template<typename derived_type, typename ord_t>
class ContainerPrintable1DBase
  : private utils::details::CrtpBase<
	ContainerPrintable1DBase<derived_type, ord_t>>{

  using this_t = ContainerPrintable1DBase<derived_type, ord_t>;
public:

  template <typename stream_t = std::ostream>
  void print(stream_t & os = std::cout,
	     char c = 'd',
	     ord_t numElementsToPrint = -1) const{
    this->underlying().printImpl(os, c, numElementsToPrint);
  }

  void printStdCout(char c = 'd', ord_t numElementsToPrint = -1) const{
    this->print(std::cout, c, numElementsToPrint);
  }

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend utils::details::CrtpBase<this_t>;

  ContainerPrintable1DBase() = default;
  ~ContainerPrintable1DBase() = default;
};//end class



template<typename derived_type, typename ord1_t, typename ord2_t = ord1_t>
class ContainerPrintable2DBase
  : private utils::details::CrtpBase<
  ContainerPrintable2DBase<derived_type, ord1_t, ord2_t>>{

  using this_t = ContainerPrintable2DBase<derived_type, ord1_t, ord2_t>;
public:

  template <typename stream_t = std::ostream>
  void print(stream_t & os = std::cout,
	     char c = 'd',
	     ord1_t niToPrint = -1,
	     ord2_t njToPrint = -1) const{
    this->underlying().printImpl(os, c, niToPrint, njToPrint);
  }

  void printStdCout(char c = 'd',
		    ord1_t niToPrint = -1,
		    ord2_t njToPrint = -1) const{
    this->print(std::cout, c, niToPrint, njToPrint);
  }

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend utils::details::CrtpBase<this_t>;

  ContainerPrintable2DBase() = default;
  ~ContainerPrintable2DBase() = default;

};//end class

}}//end namespace rompp::core
#endif
