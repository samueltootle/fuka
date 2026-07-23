#pragma once
#include <string>
#include "Configurator/config_enums.hpp"
#include "Configurator/config_utils_boost.hpp"
#include "parameter_sequence.hpp"
#include "sequence_utilities.hpp"

/**
 * \addtogroup Sequences
 * \ingroup FUKA
 * @{*/
using namespace ::Kadath::FUKA_Config;
using seq_t = Kadath::FUKA_Solvers::Parameter_sequence<BCO_PARAMS>;

namespace Kadath::FUKA_Solvers {

class ns_sequence;
std::ostream& operator<<(std::ostream&, const ns_sequence&);

class ns_sequence : public seq_t {
   protected:
    BCO_PARAMS mass_fixing_idx{BCO_PARAMS::MADM};
    double mass_fixing_val{std::nan("1")};

    BCO_PARAMS spin_fixing_idx{BCO_PARAMS::CHI};
    double spin_fixing_val{std::nan("1")};

   public:
    ns_sequence() : seq_t() {}

    /**
   * @brief Construct a new Parameter_sequence object from a string and tuple
   * objects
   *
   * @param _str Input parameter string
   * @param _t Tuple containing Config indicies for the desired parameter
   */
    ns_sequence(std::string _str, std::tuple<BCO_PARAMS> _t)
        : seq_t(_str, _t) {}

    /**
   * @brief Construct a new Parameter_sequence object from a string and a pack
   * of indices
   *
   * @param _str Input parameter string
   * @param _ts Parameter pack consisting of the desired Config indices
   */
    ns_sequence(std::string _str, BCO_PARAMS _ts) : seq_t(_str, _ts) {}

    ns_sequence(seq_t const& seq) : seq_t(seq) {}

    ns_sequence(ns_sequence const& seq) = default;

    // Getters
    double const& mass_val() const { return mass_fixing_val; }

    BCO_PARAMS const& mass_idx() const { return mass_fixing_idx; }

    double const& spin_val() const { return spin_fixing_val; }

    BCO_PARAMS const& spin_idx() const { return spin_fixing_idx; }

    // Setters
    void set_mass_val(double const& _val) { mass_fixing_val = _val; }

    void set_mass_idx(BCO_PARAMS const& _val) { mass_fixing_idx = _val; }

    void set_spin_val(double const& _val) { spin_fixing_val = _val; }

    void set_spin_idx(BCO_PARAMS const& _val) { spin_fixing_idx = _val; }

    bool is_mass_set() const { return !std::isnan(mass_fixing_val); }

    bool is_spin_set() const { return !std::isnan(spin_fixing_val); }

    std::string mass_str() const {
        auto [mass_key, mass_idx] =
            get_key_val_pair_from_val(MBCO_PARAMS, mass_fixing_idx);
        return mass_key;
    }

    std::string spin_str() const {
        auto [spin_key, spin_idx] =
            get_key_val_pair_from_val(MBCO_PARAMS, spin_fixing_idx);
        return spin_key;
    }

    BCO_PARAMS get_sequence_idx() const {
        return std::get<0>(parameter_indices);
    }

    friend std::ostream& operator<<(std::ostream&, const ns_sequence&);
};

inline bool ns_idx_is_mass_fixing(BCO_PARAMS const& idx) {
    bool is_mass_fixing{false};
    switch (idx) {
        case BCO_PARAMS::HC:
        case BCO_PARAMS::NC:
        case BCO_PARAMS::MADM:
        case BCO_PARAMS::MB:
            is_mass_fixing = true;
            break;
        default:
            break;
    }
    return is_mass_fixing;
}

inline bool ns_idx_is_spin_fixing(BCO_PARAMS const& idx) {
    bool is_spin_fixing{false};
    switch (idx) {
        case BCO_PARAMS::OMEGA:
        case BCO_PARAMS::CHI:
        case BCO_PARAMS::JADM:
            is_spin_fixing = true;
            break;
        default:
            break;
    }
    return is_spin_fixing;
}

inline void parse_tree_for_fixed(Tree const& tree,
                                 std::string const branch_name,
                                 std::string const parameter_str,
                                 BCO_PARAMS const idx,
                                 ns_sequence& seq) {
    Tree branch = read_branch(tree, branch_name);
    double _val = std::nan("1");

    extract_seq(branch, parameter_str + "_fixed", _val);

    if (!std::isnan(_val)) {
        if (ns_idx_is_mass_fixing(idx) && seq.is_mass_set()) {
            cout << std::to_string(idx) << ", "
                 << std::to_string(seq.mass_val()) << '\n';
            std::string msg{
                "Mass can only be fixed by one value. Check your config!\n"};
            throw std::runtime_error(msg.c_str());
        } else if (ns_idx_is_spin_fixing(idx) && seq.is_spin_set()) {
            std::string msg{
                "Spin can only be fixed by one value. Check your config!\n"};
            throw std::runtime_error(msg.c_str());
        } else if (ns_idx_is_mass_fixing(idx)) {
            seq.set_mass_idx(idx);
            seq.set_mass_val(_val);
        } else if (ns_idx_is_spin_fixing(idx)) {
            seq.set_spin_idx(idx);
            seq.set_spin_val(_val);
        }
    }
}

inline ns_sequence find_ns_sequence(Tree const& tree) {
    auto const& map = MBCO_PARAMS;
    std::string const branch_name = "ns";

    // Sets the Parameter Sequence up
    ns_sequence seq(find_sequence(tree, map, branch_name));

    // Need to check for fixing values
    for (const auto& [key, index] : map) {
        // Ignore resolution since this is treated separately
        if (key == "res")
            continue;
        parse_tree_for_fixed(tree, branch_name, key, index, seq);
    }
    return seq;
}

template <class seq_t>
inline bool ns_seq_is_mass_fixing(seq_t& seq) {
    auto seq_indicies = seq.get_indices();
    auto seq_idx = std::get<0>(seq_indicies);
    bool is_mass_fixing{false};
    switch (seq_idx) {
        case BCO_PARAMS::HC:
        case BCO_PARAMS::NC:
        case BCO_PARAMS::MADM:
        case BCO_PARAMS::MB:
            is_mass_fixing = true;
            break;
        default:
            break;
    }
    return is_mass_fixing;
}

template <class seq_t>
inline bool ns_seq_is_spin_fixing(seq_t& seq) {
    if (!seq.is_set())
        return false;
    auto seq_indicies = seq.get_indices();
    auto seq_idx = std::get<0>(seq_indicies);
    bool is_spin_fixing{false};
    switch (seq_idx) {
        case BCO_PARAMS::CHI:
        case BCO_PARAMS::JADM:
        case BCO_PARAMS::OMEGA:
            is_spin_fixing = true;
            break;
        default:
            break;
    }
    return is_spin_fixing;
}

inline std::ostream& operator<<(std::ostream& out, const ns_sequence& Seq) {
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
                << std::setw(20) << "step" << ": " << Seq.step_size() << '\n'
                << std::setw(20) << "# Sequences" << ": " << Seq.iterations()
                << '\n';
        } else if (Seq.is_default_set()) {
            out << std::setw(20) << "Value" << ": " << Seq.default_val()
                << '\n';
        }
        out << std::setw(20) << "indices" << ": " << indices << '\n';
    } else {
        out << "Empty Sequence\n";
    }
    s = "NS Fixing Values";
    n = ((42 - s.size()) > 0) ? 42 - s.size() : s.size() - 42;
    n /= 2;
    title = std::string(n, '*') + s + std::string(n, '*');
    out << title << std::endl;

    using namespace ::Kadath::FUKA_Config_Utils;
    auto [mass_key, mass_idx] =
        get_key_val_pair_from_val(MBCO_PARAMS, Seq.mass_idx());
    auto [spin_key, spin_idx] =
        get_key_val_pair_from_val(MBCO_PARAMS, Seq.spin_idx());
    out << std::setw(20) << mass_key << "(" << mass_idx
        << "): " << Seq.mass_val() << '\n'
        << std::setw(20) << spin_key << "(" << spin_idx
        << "): " << Seq.spin_val() << '\n'
        << std::string(42, '*') << std::endl;
    return out;
}

// Only the container values are verified - we don't modify the Config
template <class config_t>
inline void verify_ns_fixing_values(config_t& bconfig, ns_sequence& seq) {
    if (seq.is_mass_set() && seq.is_spin_set())
        return;

    if (!seq.is_mass_set()) {
        if (seq.is_set()) {
            auto idx = std::get<0>(seq.get_indices());
            if (ns_idx_is_mass_fixing(idx)) {
                seq.set_mass_idx(idx);
                seq.set_mass_val(seq.init());
            }
        }

        if (!seq.is_mass_set()) {
            if (std::isnan(bconfig.set(seq.mass_idx()))) {
                std::string msg{seq.mass_str() +
                                " value not found. Check your config.\n"};
                throw std::runtime_error(msg.c_str());
            } else {
                seq.set_mass_val(bconfig(seq.mass_idx()));
            }
        }
    }

    if (!seq.is_spin_set()) {
        if (seq.is_set()) {
            auto idx = std::get<0>(seq.get_indices());
            if (ns_idx_is_spin_fixing(idx)) {
                seq.set_spin_idx(idx);
                seq.set_spin_val(seq.init());
            }
        }

        if (!seq.is_spin_set()) {
            if (std::isnan(bconfig.set(seq.spin_idx()))) {
                std::string msg{seq.spin_str() +
                                " value not found. Check your config.\n"};
                throw std::runtime_error(msg.c_str());
            } else {
                seq.set_spin_val(bconfig(seq.spin_idx()));
            }
        }
    }
}

template <class config_t>
inline void initialize_config_from_fixing_values(config_t& bconfig,
                                                 ns_sequence const& seq) {
    bconfig.set(seq.mass_idx()) = seq.mass_val();
    bconfig.set(seq.spin_idx()) = seq.spin_val();
}

/** @}*/
}  // namespace Kadath::FUKA_Solvers
