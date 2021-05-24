#include <Rcpp.h>
#include "test.h"
using namespace Rcpp;

TestClass::TestClass(float x)
{
  this->x = x;
}

float TestClass::mult(float y)
{
  return this->x * y;
}

double test_function(TestClass c)
{
  return c.x * 2.0;
}

RCPP_EXPOSED_CLASS(TestClass);
RCPP_MODULE(test) {
  
  Rcpp::class_<TestClass>("TestClass")
  
  .constructor<float>()
  
  .field("x", &TestClass::x)
  
  .method("mult", &TestClass::mult)
  ;
  
  function("test_function", &test_function);
}

/*** R
tc <- new(TestClass, 4.0)
test_function(tc)
*/
