// version.cpp

#include "FUKA_Solvers/utilities/fuka_version.hpp"

namespace Kadath::FUKA {

std::tuple<int, int, int> parse_version(const std::string& v) {
    int major = 0, minor = 0, patch = 0;
    char dot;
    std::istringstream iss(v);
    iss >> major >> dot >> minor;
    if (iss >> dot >> patch) {
        // patch present
    }
    return {major, minor, patch};
}

bool is_older_than(const std::string& v, const std::string& target) {
    return parse_version(v) < parse_version(target);
}

std::string fuka_version() {
    return GIT_HASH;
}

std::string git_hash() {
    return GIT_HASH;
}

std::string git_describe() {
    return GIT_DESCRIBE;
}

std::string build_date() {
    return BUILD_DATE;
}

bool git_dirty() {
    return GIT_DIRTY;
}

}  // namespace Kadath::FUKA