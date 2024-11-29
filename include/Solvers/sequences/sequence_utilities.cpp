/**
 * \addtogroup Sequences
 * \ingroup FUKA
 * @{*/
using namespace Kadath::FUKA_Config;
namespace Kadath {
namespace FUKA_Solvers {

#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

template <class config_t>
config_t binary_generate_sequence_config(config_t& seqconfig,
                                         std::array<std::string, 2> bco_types,
                                         std::string outputdir) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  config_t bconfig;
  bconfig.initialize_binary(bco_types);
  bconfig.set_defaults();

  bconfig.set_filename(outputdir + "/initbin.info");

  // Copy only modified settings while retaining the remaining defaults
  // which are needed for sequence construction
  for (auto idx = 0; idx < SEQ_SETTINGS::NUM_SEQ_SETTINGS; ++idx)
    if (!std::isnan(seqconfig.seq_setting(idx)))
      bconfig.seq_setting(idx) = seqconfig.seq_setting(idx);

  // Activate all controls that are active in the seq config
  for (auto idx = 0; idx < CONTROLS::NUM_CONTROLS; ++idx)
    if (seqconfig.control(idx))
      bconfig.control(idx) = seqconfig.control(idx);

  // Activate all controls that are active in the seq config
  for (auto idx = 0; idx < STAGES::NUM_STAGES; ++idx)
    bconfig.set_stage(idx) = seqconfig.set_stage(idx);

  // Check for disabled controls that matter
  for (auto [key, idx] : MMIN_CONTROLS)
    bconfig.control(idx) = seqconfig.control(idx);

  // Copy binary parameters
  for (auto idx = 0; idx < BIN_PARAMS::NUM_BPARAMS; ++idx)
    if (!std::isnan(seqconfig.set(idx)))
      bconfig.set(idx) = seqconfig.set(idx);

  // Copy component parameters
  for (auto bco : {BCO1, BCO2}) {
    for (auto idx = 0; idx < BCO_PARAMS::NUM_BCO_PARAMS; ++idx)
      if (!std::isnan(seqconfig.set(idx, bco)))
        bconfig.set(idx, bco) = seqconfig.set(idx, bco);
  }

  // Needed for the initial imported solution for all solvers
  bconfig.control(FIXED_GOMEGA) = true;

  // make sure directory exists for outputs
  if (outputdir == "./") {
    fs::path cwd = fs::current_path();
    outputdir = cwd.string();
  }
  if (rank == 0)
    std::cout << "Binary solutions will be stored in: " << outputdir << "\n"
              << "Directory will be created if it doesn't exist.\n";
  fs::create_directory(outputdir);
  return bconfig;
}

template <class config_t>
config_t generate_sequence_config(config_t& seqconfig, std::string outputdir) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  config_t bconfig;
  bconfig.set_defaults();

  bconfig.set_filename(outputdir + "/init_co.info");

  // Copy only modified settings while retaining the remaining defaults
  // which are needed for sequence construction
  for (auto idx = 0; idx < SEQ_SETTINGS::NUM_SEQ_SETTINGS; ++idx)
    if (!std::isnan(seqconfig.seq_setting(idx)))
      bconfig.seq_setting(idx) = seqconfig.seq_setting(idx);

  // Activate all controls that are active in the seq config
  for (auto idx = 0; idx < CONTROLS::NUM_CONTROLS; ++idx)
    if (seqconfig.control(idx))
      bconfig.control(idx) = seqconfig.control(idx);

  // Check for disabled controls that matter
  for (auto [key, idx] : MMIN_CONTROLS)
    bconfig.control(idx) = seqconfig.control(idx);

  // Copy BH characteristics
  for (auto idx = 0; idx < BCO_PARAMS::NUM_BCO_PARAMS; ++idx)
    if (!std::isnan(seqconfig.set(idx)))
      bconfig.set(idx) = seqconfig.set(idx);

  // make sure BH directory exists for outputs
  if (outputdir == "./") {
    fs::path cwd = fs::current_path();
    outputdir = cwd.string();
  }
  // Activate all controls that are active in the seq config
  for (auto idx = 0; idx < STAGES::NUM_STAGES; ++idx)
    bconfig.set_stage(idx) = seqconfig.set_stage(idx);

  if (rank == 0)
    std::cout << "Isolated solutions will be stored in: " << outputdir << "\n"
              << "Directory will be created if it doesn't exist.\n";
  fs::create_directory(outputdir);
  return bconfig;
}

template <class config_t, class... bco_t>
void update_eos_parameters(config_t& seqconfig,
                           config_t& bconfig,
                           bco_t... bco) {
  for (int idx = 0; idx < EOS_PARAMS::NUM_EOS_PARAMS; ++idx)
    bconfig.set_eos(idx, bco...) = seqconfig.set_eos(idx, bco...);
}

template <class config_t, class... bco_t>
void update_diffrot_parameters(config_t& seqconfig,
                               config_t& bconfig,
                               bco_t... bco) {
  for (int idx = 0; idx < DIFFROT_PARAMS::NUM_DIFFROT_PARAMS; ++idx)
    bconfig.set_diffrot(idx, bco...) = seqconfig.set_diffrot(idx, bco...);
}

inline bool extract_seq(Tree& branch, std::string seqkey, double& storage) {
  auto [branch_name, key, val] = find_leaf(branch, seqkey);
  if (!key.empty()) {
    storage = std::stod(val);
    return true;
  }
  return false;
}

template <class... idx_t>
Parameter_sequence<idx_t...> parse_seq_tree(Tree const& tree,
                                            std::string const branch_name,
                                            std::string const parameter_str,
                                            idx_t... idx) {
  Tree branch = read_branch(tree, branch_name);
  double _val = std::nan("1"), _init = std::nan("1"), _final = std::nan("1");

  // bool vset = seqset = false;
  extract_seq(branch, parameter_str, _val);
  extract_seq(branch, parameter_str + "_init", _init);
  extract_seq(branch, parameter_str + "_final", _final);

  Parameter_sequence<idx_t...> seq(parameter_str, std::make_tuple(idx...));
  seq.set(_val, _init, _final);
  return seq;
}

template <class map_t>
uint number_of_sequences(Tree const& tree,
                         map_t const& map,
                         std::string const branch_name) {
  uint cnt{0};
  for (const auto& [key, index] : map) {
    // Ignore resolution since this is treated separately
    if (key == "res")
      continue;
    auto res = parse_seq_tree(tree, branch_name, key, index);
    if (res.is_set())
      cnt += 1;
  }
  return cnt;
}

inline uint number_of_sequences_binary(Tree const& tree) {
  uint cnt{0};
  std::string const branch_name{"binary"};
  cnt += number_of_sequences(tree, MBIN_PARAMS, branch_name);

  auto const& bco_param_map{MBCO_PARAMS};

  Tree branch = read_branch(tree, branch_name);

  auto const& node_map{MBCO};
  for (const auto& node : branch) {
    if (!node.second.empty()) {
      std::string const node_str = node.first;

      // remove suffix (e,g, bh1 -> bh, ns1 -> ns)
      auto const tstr = node_str.substr(0, 2);
      const auto& it = node_map.find(tstr);

      if (it != node_map.end())
        cnt += number_of_sequences(branch, bco_param_map, node_str);
    }
  }
  return cnt;
}

template <class map_t, class... idx_t>
decltype(auto) find_sequence(Tree const& tree,
                             map_t const& map,
                             std::string const branch_name,
                             idx_t... idx) {
  for (const auto& [key, index] : map) {
    // Ignore resolution since this is treated separately
    if (key == "res")
      continue;
    auto res = parse_seq_tree(tree, branch_name, key, index, idx...);
    if (res.is_set())
      return res;
  }
  return parse_seq_tree(tree, branch_name, "", (*map.begin()).second, idx...);
}

inline decltype(auto) find_sequence_binary(Tree const& tree) {
  std::string const branch_name{"binary"};

  auto const& bco_param_map{MBCO_PARAMS};

  Tree branch = read_branch(tree, branch_name);

  auto const& node_map{MBCO};

  std::array<NODES, 2> bcos{NODES::BCO1, NODES::BCO2};
  for (const auto& node : branch) {
    if (!node.second.empty()) {
      std::string const node_str = node.first;

      // remove suffix (e,g, bh1 -> bh, ns1 -> ns)
      auto const tstr = node_str.substr(0, 2);

      const auto& it = node_map.find(tstr);

      if (it != node_map.end()) {
        // save suffix {1, 2}
        int const tstr_suffix = std::atoi(&node_str.back());
        if (tstr_suffix > 2 || tstr_suffix < 1) {
          std::cerr << node_str << " not recognized in Config\n";
          std::_Exit(EXIT_FAILURE);
        }
        int const idx = tstr_suffix - 1;
        auto res = find_sequence(branch, MBCO_PARAMS, node_str, bcos[idx]);
        return res;
      }
    }
  }
  return parse_seq_tree(tree, branch_name, "", (*bco_param_map.begin()).second,
                        BCO1);
}

template <class TupType, size_t... I>
std::ostream& tuple_print(std::ostream& os,
                          const TupType& _tup,
                          std::index_sequence<I...>) {
  os << "(";
  (..., (os << (I == 0 ? "" : ", ") << std::get<I>(_tup)));
  os << ")";
  return os;
}

template <class... Ts>
std::ostream& operator<<(std::ostream& out,
                         const Parameter_sequence<Ts...>& Seq) {
  auto indices = Seq.get_indices();
  std::string s = Seq.str() + " sequence";
  int n = ((42 - s.size()) > 0) ? 42 - s.size() : s.size() - 42;
  n /= 2;
  std::string title = std::string(n, '*') + s + std::string(n, '*');
  out << title << std::endl;
  if (Seq.is_set() || Seq.is_default_set()) {
    if (Seq.is_set()) {
      out << std::setw(20) << "initial" << ": " << Seq.init() << '\n'
          << std::setw(20) << "final" << ": " << Seq.final() << '\n'
          << std::setw(20) << "step" << ": " << Seq.step_size() << '\n';
      if (Seq.str() != "res")
        out << std::setw(20) << "# Sequences" << ": " << Seq.iterations()
            << '\n';
    } else if (Seq.is_default_set()) {
      out << std::setw(20) << "Value" << ": " << Seq.default_val() << '\n';
    }
    out << std::setw(20) << "indices" << ": " << indices << '\n';
  } else {
    out << "Empty Sequence\n";
  }
  return out;
}

template <class... T>
std::ostream& operator<<(std::ostream& os, const std::tuple<T...>& _tup) {
  return tuple_print(os, _tup, std::make_index_sequence<sizeof...(T)>());
}

template <class config_t, class Res_t>
void verify_resolution_sequence(config_t& bconfig, Res_t& resolution) {
  if (resolution.is_set())
    return;

  // Determine lowest resolution
  auto init_res = (std::isnan(bconfig.seq_setting(SEQ_SETTINGS::INIT_RES)))
                      ? 9
                      : bconfig.seq_setting(SEQ_SETTINGS::INIT_RES);

  // Determine highest resolution
  auto final_res =
      (resolution.is_default_set()) ? resolution.default_val() : init_res;
  if (final_res < init_res)
    std::swap(init_res, final_res);

  resolution.set(final_res, init_res, final_res);
}

template <class config_t, class seq_t>
void update_filename_from_mass_fixing(config_t& bconfig,
                                      seq_t& seq,
                                      std::stringstream& ss) {
  auto idx{seq->mass_idx()};
  auto [seq_key, tidx] = get_key_val_pair_from_val(MBCO_PARAMS, idx);

  switch (idx) {
    case BCO_PARAMS::HC:
    case BCO_PARAMS::NC:
    case BCO_PARAMS::MADM:
    case BCO_PARAMS::MB:
      ss << seq_key << "." << bconfig(idx) << ".";
      break;
    default:
      std::string msg{
          "Sequence initialized, but not implemented for Mass index = " +
          std::to_string(int(idx))};
      throw std::runtime_error(msg.c_str());
      break;
  }
}

template <class config_t, class seq_t>
void update_filename_from_spin_fixing(config_t& bconfig,
                                      seq_t& seq,
                                      std::stringstream& ss) {
  auto idx{seq->spin_idx()};
  auto [seq_key, tidx] = get_key_val_pair_from_val(MBCO_PARAMS, idx);

  switch (idx) {
    case BCO_PARAMS::OMEGA:
    case BCO_PARAMS::JADM:
    case BCO_PARAMS::CHI:
      ss << seq_key << "." << bconfig(idx) << ".";
      break;
    default:
      std::string msg{
          "Sequence initialized, but not implemented for Mass index = " +
          std::to_string(int(idx))};
      throw std::runtime_error(msg.c_str());
      break;
      break;
  }
}

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath
