#pragma once
#include <string>

#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

namespace Kadath::FUKA_Solvers {

    /**
     * @brief Get the global path to saved compact object solutions (COs)
     *
     * @return std::string
     */
    inline std::string get_cos_path() {
        const std::string home_kadath{std::getenv("HOME_KADATH")};
        std::string central_abs{home_kadath + "/COs/"};
        return central_abs;
    }
}