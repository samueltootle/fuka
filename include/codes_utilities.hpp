//
// Created by sauliac on 29/05/2020.
//
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "headcpp.hpp"

#ifndef __KADATH_CODES_UTILITY_HPP_
#define __KADATH_CODES_UTILITY_HPP_

//! Macro to declare a pointer data member with associated trivial accessors.
#define ptr_data_member(type,identifier,smart_ptr_type) \
protected:\
    std:: smart_ptr_type##_ptr<type> identifier;\
public:\
    std:: smart_ptr_type##_ptr<type> const & get_##identifier() const {return identifier;}\
    std:: smart_ptr_type##_ptr<type> & get_##identifier() {return identifier;}

//! Macro to declare internal variable with a read-only trivial accessor.
#define internal_variable(type,identifier) \
protected:\
    type identifier;\
public:\
    optimal_access_type<type> get_##identifier () const {return identifier;}



/****************************************************************************************
 * A simple argument parser.
 ****************************************************************************************/

template<typename T> struct is_admissible_option_type {
    static constexpr bool value = false;
};
template<> struct is_admissible_option_type<bool> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<int> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<long> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<long long> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<double> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<std::string> {static constexpr bool value = true;};
template<> struct is_admissible_option_type<char*> {static constexpr bool value = true;};

template<typename T> struct option_type_name {static std::string get() {return typeid(T).name();}};
template<> struct option_type_name<bool> {static std::string get() {return "bool";}};
template<> struct option_type_name<int> {static std::string get() {return "int";}};
template<> struct option_type_name<long> {static std::string get() {return "long";}};
template<> struct option_type_name<long long> {static std::string get() {return "long long";}};
template<> struct option_type_name<double> {static std::string get() {return "double";}};
template<> struct option_type_name<std::string> {static std::string get() {return "string";}};
template<> struct option_type_name<char*> {static std::string get() {return "char*";}};


template<typename T> inline T from_string(std::string const & s) {
    std::stringstream sss{s};
    T t_from_s{};
    sss >> t_from_s;
    return t_from_s;
}
template<> inline bool from_string<bool>(std::string const & s) {
    if( s == "0" || s == "false" || s == "False" || s == "FALSE" ||
        s == "off" || s == "Off" || s == "OFF" || s == "disabled" ||
        s == "Disabled" || s == "DISABLED") {
        return false;
    }
    else if( s == "1" || s == "true" || s == "True" || s == "TRUE" ||
        s == "on" || s == "On" || s == "ON" || s == "enabled" ||
        s == "Enabled" || s == "ENABLED") {
        return true;
    }
    else throw std::runtime_error{"bad converion from string too bool"};
}
template<> inline std::string from_string<std::string>(std::string const & s) {return s;}
template<> inline char const * from_string<char const *>(std::string const & s) {return s.c_str();}

std::vector<std::string> inline part_string(std::string str,std::size_t part_size) {
    std::vector<std::string> parts{};
    while(!str.empty()) {
        auto pos = str.find_first_of(' ',part_size);
//        auto pos = std::find(str.cbegin()+part_size,str.cend(),' ');
        parts.push_back(str.substr(0,pos));
        str.erase(0,pos);
        pos = str.find_first_not_of(' ');
        str.erase(0,pos);
    }
    return std::move(parts);
}

struct Option_base {
    static constexpr std::size_t tab_size{4};
    static constexpr std::size_t max_key_size{15};
    static constexpr std::size_t max_type_name_size{15};
    static constexpr std::size_t description_tab_size{3*tab_size + max_key_size + max_type_name_size};
    static constexpr std::size_t description_column_size {100};
    std::string key;
    std::string type_name;
    std::string description;
    Option_base() = default;
    Option_base(std::string _key, std::string _type_name,std::string _description) :
        key{std::move(_key)}, type_name{std::move(_type_name)}, description{std::move(_description)} {}
    virtual ~Option_base() = default;
    virtual Option_base & set(std::string const & input) = 0;
    virtual void display(std::ostream &os) const {
        auto const parted_descr = part_string(description,description_column_size);
        os << "    " << std::setw(max_key_size) << key;
        if(key=="-h") os << "    " << std::setw(max_type_name_size) << "N/A" << "    ";
        else os << "    " << std::setw(max_type_name_size) << type_name << "    ";
        bool first_line {true};
        for(auto const & line : parted_descr) {
            if(!first_line) {
                os << std::setw(description_tab_size ) << " " << line << '\n';
            }
            else {
                os << line << '\n';
                first_line = false;
            }
        }
    }
};
std::ostream & operator<<(std::ostream & os,Option_base const & option) {option.display(os); return os;}

template<typename T> struct Option : Option_base {
    T value;
    T default_value;
    CXX_17_ATTRIBUTES(maybe_unused) Option(std::string const & key,std::string const & description,T _default_value,std::string const & input)
        : Option_base{key,option_type_name<T>::get(),description}, value{}, default_value{_default_value}
    {
        this->set(input);
    }
    Option() = default;
    Option & set(std::string const & input) override {
        if(input.empty()) value = default_value;
        else value = from_string<T>(input);
        return *this;
    }
    void display(std::ostream & os) const override {
        this->Option_base::display(os);
        os << std::setw(description_tab_size) << ' ';
        os << "Default value : " << default_value << std::endl;
    }
};

template<> struct Option<bool> : Option_base {
    bool value{};
    bool default_value{};
    CXX_17_ATTRIBUTES(maybe_unused) Option(std::string const & key,std::string const & description,bool _default_value,
            std::string const & input = "") :
        Option_base{key,"bool",description}, value{}, default_value{_default_value}
    {}
    Option() = default;
    Option & set(std::string const & s) override {
        if(s.empty()) value = default_value;
        else value = from_string<bool>(s);
        return *this;}
    void display(std::ostream & os) const override {
        this->Option_base::display(os);
        if(key != "-h") {
            os << std::setw(description_tab_size) << ' ';
            os << "Default value : " << (default_value ? "true" : "false") << std::endl;
        }
    }
};


class Arguments_parser{
public:
private:
    std::string executable;
    std::vector<std::string> command_line;
    std::map<std::string,std::unique_ptr<Option_base>> option_list;

public:
    Arguments_parser (int &argc, char **argv,std::string _executable="") :
        executable{std::move(_executable)}, command_line{}, option_list{} {
        if(executable.empty()) executable = argv[0];
        for (int i=1; i < argc; ++i)
            this->command_line.emplace_back(argv[i]);
    }

    CXX_17_ATTRIBUTES(nodiscard) std::string const & get_executable() const {return executable;}
    Arguments_parser & set_executable(std::string const & new_executable) {executable = new_executable; return *this;}

    template<typename T>
    CXX_17_ATTRIBUTES(maybe_unused) void reference_option(std::string const &key,std::string const &description,T default_value,
                                                std::string const & input= "") {
        std::unique_ptr<Option_base> opt{new Option<T>{key,description,default_value,input}};
        auto pos = option_list.find(key);
	    if(pos == option_list.end()) option_list.emplace(key,std::move(opt));
    }
    void reference_option(std::string const &key,std::string const &description,bool value) {
        std::unique_ptr<Option_base> opt{new Option<bool>{key,description,value}};
        auto pos = option_list.find(key);
	    if(pos == option_list.end()) option_list.emplace(key,std::move(opt));
    }

    bool find_option(const std::string &option,std::string const &description) {
        auto itr = std::find(command_line.cbegin(),command_line.cend(),option);
        bool const found{itr != command_line.end()};
        reference_option(option,description,found);
        return found;
    }

    template<typename OT>
    std::pair<typename std::enable_if<
            is_admissible_option_type<OT>::value,
            OT
    >::type,bool>
    get_option_value(const std::string &option,std::string const &description,OT default_value)  {
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->command_line.begin(), this->command_line.end(), option);
        reference_option(option,description,default_value,(itr != command_line.end() ? *itr : "" ));
        if (itr != this->command_line.end() && ++itr != this->command_line.end()){
            return {from_string<OT>(*itr),true};
        }
        else return {default_value,false};
    }

    template<typename OT>
    std::pair<typename std::enable_if<
                is_admissible_option_type<OT>::value,
                OT
            >::type,bool>
    get_option_value(const std::string &option) const {
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->command_line.begin(), this->command_line.end(), option);
        if (itr != this->command_line.end() && ++itr != this->command_line.end()){
            return {from_string<OT>(*itr),true};
        }
        else return {OT{},false};
    }

    virtual void display(std::ostream & os) const {
        os << executable << " executable, usage : " << executable << " [options] " << '\n';
        os << "Available options descriptions : " << std::endl;
        os << "    " << std::setw(Option_base::max_key_size) << "   Option Key  "
           << "    " << std::setw(Option_base::max_type_name_size) << "  Option Type  "
           << "    " << std::setw(Option_base::description_column_size/2 + 6) << "Description" << std::endl;
        for(auto const & option : option_list) {
            option.second->display(os);

        }
    }
};
std::ostream& operator<<(std::ostream& os,Arguments_parser const & arguments_parser) {
    arguments_parser.display(os);
    return os;
}



#endif //__KADATH_CODES_UTILITY_HPP_
