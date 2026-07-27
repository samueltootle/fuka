// version.cpp

#include "FUKA_Solvers/utilities/fuka_version.hpp"

namespace Kadath::FUKA {

std::string git_hash() {
    return GIT_HASH;
}

std::string git_describe() {
    return GIT_DESCRIBE;
}

bool git_dirty() {
    return GIT_DIRTY;
}

}  // namespace Kadath::FUKA