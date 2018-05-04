#ifndef EVECTOR_HPP
#define EVECTOR_HPP

#include <ostream>
#include "epetramock_config.h"
#include <vector>


namespace epetramock
{

class evector
{
public:
	evector(){}
	~evector(){}

	void print();
private:
	int mySize_;
	std::vector<double> data_;	
};

/** \brief . */
// std::string deps();

/* \brief Simple hello world class.
 */
// class HelloWorld {
// public:
//   /** \brief. */
//   HelloWorld();
//   /** \brief . */
//   void printHelloWorld(std::ostream &out) const;
//   * \brief Deprecated. 
//   SIMPLECXX_DEPRECATED int someOldFunc() const; 
//   /** \brief Deprecated. */
//   SIMPLECXX_DEPRECATED_MSG("Just don't call this function at all please!")
//   int someOldFunc2() const; 
// };


} // namespace
#endif
