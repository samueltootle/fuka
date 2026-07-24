#pragma once
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <variant>
#include "name_tools.hpp"
namespace pt = boost::property_tree;
using Tree = pt::ptree;

/**
 * \addtogroup Configurator
 * @{*/
namespace Kadath {
namespace FUKA_Config {
// using namespace Kadath::FUKA_Config_Utils;

class configurator_base {
   protected:
    std::string filename{};
    std::string outputdir{"./"};
    Tree tree;

   public:
    configurator_base() = default;
    configurator_base(std::string ifile) : filename(ifile) {};

    /**
   * configurator_base::config_filename_abs
   *
   * returns the config filename based with absolute path.  This reduced
   * a lot of repeat code.
   *
   * @return string with <outputdir/filename.info>
   */
    const std::string config_filename_abs() const;

    /**
   * configurator_base::set_filename
   *
   * Sets the config filename.  If the input string contains a path '/'
   * then the output directory is also updated.
   * Will also automatically append .info if missing.
   *
   * @param[input] filename: string with filename and possibly the output dir
   */
    void set_filename(std::string fname);

    /**
   * configurator_base::set_outputdir
   *
   * Sets the config output directory.  Appends '/' if not found
   *
   * @param[input] dir: string with the output dir
   */
    void set_outputdir(std::string dir);

    /**
   * configurator_base::config_filename()
   * @return filename: returns constant filename
   */
    const std::string config_filename() const { return filename; }

    void read_config();

    Tree const& get_config_tree() const { return tree; }
};

/**
 * @}*/
}  // namespace FUKA_Config
}  // namespace Kadath