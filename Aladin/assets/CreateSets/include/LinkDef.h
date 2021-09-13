#include <vector>
#include <map>
#ifdef __ROOTCLING__
#pragma link C++ class std::map<double, std::vector<double>>+;
#pragma link C++ class std::map<double, double>+;
#endif