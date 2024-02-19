#pragma once
#include<string>
#include "Configurator/config_binary.hpp"
#include "Configurator/configurator_boost.hpp"
namespace Kadath::FUKA_Solvers {
template<typename config_t, typename... bco_idx>
/**
 * @brief Extract EOS name from the Config object
 * 
 * @tparam bco_idx optional index type for binary Configs
 * @param bco optional index for binary Configs
 * @return std::string 
 */
template<typename... bco_idx>
std::string extract_eos_name(config_t& bconfig, bco_idx... bco) {
const std::string eos_file_abs = bconfig.template eos<std::string>(EOSFILE, bco...);
const std::string eos_file = extract_filename(eos_file_abs);
return eos_file.substr(0, eos_file.find("."));
}
}