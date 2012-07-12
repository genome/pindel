#ifndef FN_PARAMETERS_H
#define FN_PARAMETERS_H

#include <string>
#include <vector>

#include "parameter.h"

void printHelp(const std::vector<Parameter *>& parameters );
void readParameters(int argc, char *argv[], std::vector<Parameter *>& parameters);
void defineParameters( std::vector<Parameter *>& parameters );
unsigned int findParameter(const std::string& name, const std::vector<Parameter *>& parameters);
bool checkParameters(const std::vector<Parameter *>& parameters);

#endif
