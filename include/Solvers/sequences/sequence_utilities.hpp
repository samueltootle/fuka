#pragma once
#include "Configurator/config_binary.hpp"
#include "parameter_sequence.hpp"
#include<mpi.h>
#include<string>
#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
using namespace std::experimental;
#else
#include <filesystem>
#endif
#include<algorithm>

/**
 * \addtogroup Sequences
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

/**
 * @brief Update Config EOS parameters based on sequence Config
 * 
 * @tparam config_t Sequence Config type
 * @tparam bco_t Optional pack for compact object to read from
 * @param seqconfig Sequence Config
 * @param bconfig Full Config
 * @param bco Optional index for the relevant CO to update
 */
template<class config_t, class... bco_t>
void update_eos_parameters (config_t & seqconfig, config_t& bconfig, bco_t... bco);

/**
 * @brief Generate full Config file from a minimal sequence Config
 * 
 * @tparam config_t Config type
 * @param seqconfig Sequence Config
 * @param outputdir Output location
 * @return config_t Output directory
 */
template<class config_t>
config_t generate_sequence_config (config_t & seqconfig, std::string outputdir);

/**
 * @brief Generate full binary Config from a minimal sequence Config
 * 
 * @tparam config_t Config type
 * @param seqconfig Sequence Config
 * @param bco_types Object types <"ns" or "bh">
 * @param outputdir Output location
 * @return config_t Output directory
 */
template<class config_t>
config_t binary_generate_sequence_config (
    config_t & seqconfig, std::array<std::string, 2> bco_types, std::string outputdir);

/**
 * @brief Attempt to extract a sequence parameter from a BOOST tree
 * 
 * @param branch Branch to look for a sequence
 * @param seqkey Key of the branch leaf to look for
 * @param storage location to store the value of the result
 * @return true 
 * @return false 
 */
inline bool extract_seq(Tree& branch, std::string seqkey, double& storage);

/**
 * @brief Build a parameter sequence based on the input parameters parsed from
 * a tree
 * 
 * @tparam idx_t Possibly stored index enum types
 * @param tree  Boost tree to search
 * @param branch_name Branch to search for sequences
 * @param parameter_str Parameter to search for
 * @param idx Pack of indicies to be stored in the parameter_sequence
 * @return Parameter_sequence Derived sequence
 */
template<class... idx_t>
Parameter_sequence<idx_t...> parse_seq_tree(Tree const & tree, 
    std::string const branch_name, std::string const parameter_str, idx_t... idx);

/**
 * @brief Non-recursive search for number of detected sequences based on the given inputs - needed to verify
 * only one sequence is present.
 * 
 * @tparam map_t std::map<string, ENUM> where ENUM
 * @param tree  Boost tree to search
 * @param map Map containing potential inputs
 * @param branch_name Node to read from tree (binary, ns, etc) 
 * @return uint Number of sequences found
 */
template<class map_t>
uint number_of_sequences(Tree const & tree, map_t const & map, std::string const branch_name);

/**
 * @brief Static recursive search for number of detected sequences 
 * based on the given inputs - needed to verify only one sequence is present.
 * 
 * @param tree  Boost tree to search
 * @return uint Number of sequences found
 */
inline uint number_of_sequences_binary(Tree const & tree);

/**
 * @brief Search for a sequence in a tree based on the inputs
 * 
 * @tparam map_t std::map<string, ENUM> where ENUM
 * @tparam idx_t pack of potential ENUM types
 * @param tree  Boost tree to search
 * @param map Map containing potential inputs
 * @param branch_name Node to read from tree (binary, ns, etc)
 * @return Parameter_sequence Empty or discovered Parameter sequence
 */
template<class map_t, class... idx_t>
decltype(auto) find_sequence(Tree const & tree, map_t const & map, 
  std::string const branch_name, idx_t... idx);

/**
 * @brief Recursive search for a sequence in a tree based on the inputs
 * 
 * @param tree  Boost tree to search
 * @return uint Number of sequences found
 */
inline decltype(auto) find_sequence_binary(Tree const & tree);

/**
 * @brief Repurposed overloads from stackoverflow to print tuples nicely
 * 
 * @tparam TupType tuple type
 * @tparam I Index sequence to unpack tuple
 * @param os output stream
 * @param _tup input tuple
 * @return std::ostream& 
 */
template<class TupType, size_t... I>
std::ostream& tuple_print(std::ostream& os,
                          const TupType& _tup, std::index_sequence<I...>);

/**
 * ostream operator<<
 * Overloaded operator<< to handle printing Parameter Sequences to standard output
 *
 * @tparam Ts Enum types
 * @param out: output stream
 * @param seq: Parameter sequence
 * @return std::ostream& 
 */
template<class... Ts>
std::ostream& operator<<(std::ostream& out, const Parameter_sequence<Ts...>& Seq);

/**
 * ostream operator<<
 * Overloaded operator<< to handle printing Resolution Sequences to standard output
 *
 * @tparam Ts Enum types
 * @param out: output stream
 * @param seq: Resolution sequence
 * @return std::ostream& 
 */
template<class... T>
std::ostream& operator<< (std::ostream& os, const std::tuple<T...>& _tup);

/**
 * @brief Verify resolution sequence for a new ID sequence. Set defaults
 * if the default_value is larger than the default of 9pts.  No action
 * if a valid sequence is passed
 * 
 * @tparam config_t Config type
 * @tparam Res_t Resolution Sequence type
 * @param bconfig Configurator object
 * @param resolution Resolution Sequence object
 */
template<class config_t, class Res_t>
void verify_resolution_sequence(config_t& bconfig, Res_t& resolution);

/**
 * @brief Modify the output filename based on the initialized sequence
 * FIXME: This is currently configured for NS initial data.  Needs to be
 * more generic
 * 
 * @tparam config_t Configurator type
 * @tparam seq_t Sequence type
 * @param bconfig Config
 * @param seq Sequence
 * @param ss output filename stringstream
 */
template<class config_t, class seq_t>
void update_filename_from_mass_fixing(config_t& bconfig, seq_t & seq, std::stringstream& ss);

/**
 * @brief Modify the output filename based on the initialized sequence
 * FIXME: This is currently configured for NS initial data.  Needs to be
 * more generic
 * 
 * @tparam config_t Configurator type
 * @tparam seq_t Sequence type
 * @param bconfig Config
 * @param seq Sequence
 * @param ss output filename stringstream
 */
template<class config_t, class seq_t>
void update_filename_from_spin_fixing(config_t& bconfig, seq_t& seq, std::stringstream& ss);

/**
 * @brief Determine if a sequence is based on fixing the mass
 * 
 * @tparam seq_t Sequence type
 * @param seq sequence
 * @return bool
 */
template<class seq_t>
bool seq_is_mass_fixing(seq_t& seq);
/** @}*/
}}
#include "sequence_utilities.cpp"