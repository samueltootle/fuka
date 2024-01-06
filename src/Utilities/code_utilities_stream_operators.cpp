#include <codes_utilities.hpp>

std::ostream & operator<<(std::ostream & os,Option_base const & option) {option.display(os); return os;}

std::ostream& operator<<(std::ostream& os,Arguments_parser const & arguments_parser) {
  arguments_parser.display(os);
  return os;
}
