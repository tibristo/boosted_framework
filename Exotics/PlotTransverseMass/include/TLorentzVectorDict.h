#include <vector>
#include "TLorentzVector.h"

#ifdef __CINT__
#pragma link C++ class std::vector< TLorentzVector >;
#pragma link C++ class std::vector< float >;
#pragma link C++ class std::vector< std::vector< int > >;
#endif

template class std::vector< TLorentzVector >;
template class std::vector< float >;
template class std::vector< std::vector< int > >;

