#pragma once

#include <exception>
#include <iostream>

class InputException: public std::exception {
public:
  InputException(std::string msg) : msg_(msg){};
  ~InputException() throw() {};
  const char* what() const throw() { return msg_.c_str();}

private:
  std::string msg_;
};