#include <Configurator/configurator_base.hpp>

namespace Kadath {
namespace FUKA_Config {
const std::string configurator_base::config_filename_abs() const {
    int idx = filename.rfind(".");
    return std::string{outputdir + filename};
}

void configurator_base::set_filename(std::string fname) {
    if (fname.find('/') == std::string::npos)
        filename = fname;
    else {
        filename = Kadath::extract_filename(fname);
        outputdir = Kadath::extract_path(fname);
    }
    if (filename.find(".info") == std::string::npos)
        filename += ".info";
    return;
}

void configurator_base::set_outputdir(std::string dir) {
    if (dir.back() != '/')
        dir.push_back('/');
    this->outputdir = std::move(dir);
}

void configurator_base::read_config() {
    pt::read_info(this->config_filename_abs(), tree);
};
}  // namespace FUKA_Config
}  // namespace Kadath