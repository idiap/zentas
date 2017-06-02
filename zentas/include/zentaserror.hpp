#ifndef ZENTASERROR_HPP
#define ZENTASERROR_HPP

#include <stdexcept>

namespace zentas{

class zentas_error : public std::runtime_error{
  public:
    zentas_error(const std::string& what_arg);
};

void zentas_warning(const std::string & warning);

}


#endif
