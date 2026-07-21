/*
 * Copyright 2022
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 * L. Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

// build a shared library for interacting with the Configurator files in python
#include "Configurator/pyconfigurator.hpp"

using namespace Kadath;
using namespace Kadath::FUKA_pyTools;
using namespace Kadath::FUKA_Config;

BOOST_PYTHON_MODULE(_pyConfigurator) {
    // Add reader for binary config files
    Kadath::FUKA_pyTools::constructPythonConfigurator<
        bin_Configurator_reader_t>("Binary_Configurator");

    using bh_config_t = kadath_config_boost<BCO_BH_INFO>;
    using bh_Configurator_reader_t = Configurator_reader_t<bh_config_t>;
    // Add reader for Generic compact object config files
    Kadath::FUKA_pyTools::constructPythonConfigurator<bh_Configurator_reader_t>(
        ("BH_Configurator"));

    // Add reader for Generic compact object config files
    Kadath::FUKA_pyTools::constructPythonConfigurator<ns_Configurator_reader_t>(
        ("NS_Configurator"));
}