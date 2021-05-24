#include <Rcpp.h>
using namespace Rcpp;

#ifndef _TEST_H
#define _TEST_H

class TestClass
{
public:
  TestClass(float x);
  
  
  float mult(float y);
  
  float x;
  
};

double test_function(TestClass c);

#endif