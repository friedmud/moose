#include "ExceptionSteady.h"

template<>
InputParameters validParams<ExceptionSteady>()
{
  return validParams<Steady>();
}

ExceptionSteady::ExceptionSteady(const std::string & name, InputParameters parameters) :
    Steady(name, parameters)
{
}

ExceptionSteady::~ExceptionSteady()
{
}

void
ExceptionSteady::execute()
{
}
