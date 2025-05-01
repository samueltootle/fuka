/*
 * This file is part of the KADATH library.
 * Author: Samuel Tootle
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

#pragma once
#include <algorithm>
#include "configurator_boost.hpp"
#include "name_tools.hpp"

/**
 * Default constructor for a given Parameter Container
 *
 * @tparam ParamC type of parameter container (e.g. BH_INFO, BIN_INFO)
 */
template <typename ParamC>
kadath_config_boost<ParamC>::kadath_config_boost() : configurator_base() {
  fields.fill(false);
  stages.fill(false);
  controls.fill(false);
  seq_settings.fill(std::nan("1"));
  set_seq_defaults();
}

template <typename ParamC>
kadath_config_boost<ParamC>::kadath_config_boost(std::string ifile)
    : configurator_base(ifile), container{} {

  this->set_filename(ifile);
  this->open_config();
  set_seq_defaults();
}

template <typename ParamC>
int kadath_config_boost<ParamC>::open_config() {
  read_config();
  int status = -1;
  auto key = MBCO.find(container.get_type());
  if (container.get_type() == "binary" || key != MBCO.end()) {
    branch_in = read_branch(tree, container.get_type());
  } else {
    for (auto& e : MBCO) {
      branch_in = read_branch(tree, e.first);
      if (!branch_in.empty()) {
        status = e.second;
        break;
      }
    }
  }
  if (branch_in.empty()) {
    std::cout << "No Node found. \n";
    return status;
  }
  container.read_params(branch_in);

  if (tree.find("fields") != tree.not_found())
    read_keys(MBCO_FIELDS, fields, read_branch(tree, "fields"));
  if (tree.find("stages") != tree.not_found())
    read_keys(MSTAGE, stages, read_branch(tree, "stages"));
  if (tree.find("sequence_controls") != tree.not_found())
    read_keys(MCONTROLS, controls, read_branch(tree, "sequence_controls"));
  if (tree.find("sequence_settings") != tree.not_found())
    read_keys(MSEQ_SETTINGS, seq_settings,
              read_branch(tree, "sequence_settings"));
  return status;
}

template <typename ParamC>
void kadath_config_boost<ParamC>::write_config(std::string ofile) {

  Tree new_tree;
  Tree branch = container.return_branch();

  std::string s{"initial"};
  if (tree.find(s) != tree.not_found()) {
    if (this->control(UPDATE_INIT))
      new_tree.put_child(s, branch_in);
    else
      new_tree.push_back(std::make_pair(s, tree.get_child(s)));
  } else {
    if (this->control(UPDATE_INIT))
      new_tree.put_child(s, branch);
  }

  // get branch from parameter container
  s = container.get_type();
  for (auto& key : branch) {
    if (key.first == s) {
      new_tree.push_back(std::make_pair(s, key.second));
    } else {
      std::string news = s + "." + key.first;
      if (key.second.empty())
        new_tree.put(news.c_str(), key.second.data());
      else
        add_branch_data(new_tree, key.second, news);
    }
  }

  // add branches for fields, stages, and controls
  s = "fields";
  new_tree.push_back(
      std::make_pair(s, build_branch<Tree>(MBCO_FIELDS, fields)));
  s = "stages";
  auto stage_map = append_map(MSTAGE, container.get_stage_map(), stages);
  new_tree.push_back(
      std::make_pair(s, build_branch<Tree>(stage_map, stages, true)));
  s = "sequence_controls";
  new_tree.push_back(
      std::make_pair(s, build_branch<Tree>(MCONTROLS, controls, true)));
  s = "sequence_settings";
  new_tree.push_back(
      std::make_pair(s, build_branch<Tree>(MSEQ_SETTINGS, seq_settings)));

  if (ofile != "null") {
    //this will update outputdir if ofile includes a path
    set_filename(ofile);
  }

  s = config_filename_abs();
  pt::write_info(s, new_tree);
}

template <typename ParamC>
void kadath_config_boost<ParamC>::write_minimal_config(std::string ofile) {

  Tree new_tree;
  Tree branch = container.return_branch();
  std::string s{"initial"};

  // get branch from parameter container
  s = container.get_type();
  for (auto& key : branch) {
    if (key.first == s) {
      new_tree.push_back(std::make_pair(s, key.second));
    } else {
      std::string news = s + "." + key.first;
      if (key.second.empty())
        new_tree.put(news.c_str(), key.second.data());
      else
        add_branch_data(new_tree, key.second, news);
    }
  }

  // add branches for fields, stages, and controls
  s = "stages";
  auto stage_map = append_map(MSTAGE, container.get_stage_map(), stages);
  new_tree.push_back(
      std::make_pair(s, build_branch<Tree>(stage_map, stages, true)));
  s = "sequence_settings";
  new_tree.push_back(
      std::make_pair(s, build_branch<Tree>(MSEQ_SETTINGS, seq_settings)));
  s = "sequence_controls";
  auto control_map = append_map(MCONTROLS, MMIN_CONTROLS, controls);
  new_tree.push_back(
      std::make_pair(s, build_branch<Tree>(control_map, controls, true)));

  if (ofile != "null") {
    //this will update outputdir if ofile includes a path
    set_filename(ofile);
  }

  s = config_filename_abs();
  pt::write_info(s, new_tree);
}

/**
 * ostream operator<<
 * Overloaded operator<< to handle printing configuration file to standard output
 *
 * @tparam T type of parameter container (e.g. BH_INFO, BIN_INFO)
 * @param[input] out: output stream
 * @param[input] config: Configurator config to print
 */
template <typename T>
std::ostream& operator<<(std::ostream& out,
                         const kadath_config_boost<T>& config) {
  out << config.container;
  return out;
}

/**
 * kadath_config_boost::space_filename
 *
 * returns the space filename based on the config filename.  This reduced
 * a lot of repeat code.
 * 
 * @return string with <outputdir/filename.dat>
 */
template <typename ParamC>
const std::string kadath_config_boost<ParamC>::space_filename() const {
  int idx = filename.rfind(".");
  return std::string{outputdir + filename.substr(0, idx) + ".dat"};
}

template <typename ParamC>
inline void kadath_config_boost<ParamC>::set_seq_defaults() {
  if (std::isnan(seq_settings[PREC]))
    seq_settings[PREC] = 1e-8;
  if (std::isnan(seq_settings[MAX_ITER]))
    seq_settings[MAX_ITER] = 15;
}
