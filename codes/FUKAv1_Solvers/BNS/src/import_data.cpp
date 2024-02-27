#include "kadath.hpp"
#include "Configurator/config_binary.hpp"

#include <sstream>
#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif
using namespace Kadath::FUKA_Config;

int main(int argc, char** argv) {

	// read the input NS'
  if(argc < 4) {
    std::cerr << "Missing input config files \n i.e. ./import_data binary.info ns1.info ns2.info";
    std::_Exit(EXIT_FAILURE);
  }

	std::stringstream ss;

  if(argc < 5)
    ss << "import";
  else
    ss << argv[4];

  std::cout << "Importing data into this binary:" << std::endl;

  // read in binary configuration
  std::string bin_fn = argv[1];
  kadath_config_boost<BIN_INFO> bconfig(bin_fn);

  std::cout << bconfig << std::endl;

  std::cout << "Importing data from these stars:" << std::endl;

  // read in single NS configurations
  std::string ns1_fn{argv[2]};
  std::string ns2_fn{argv[3]};

  kadath_config_boost<BCO_NS_INFO> ns1_config(ns1_fn);
  kadath_config_boost<BCO_NS_INFO> ns2_config(ns2_fn);

  std::cout << ns1_config << std::endl;
  std::cout << ns2_config << std::endl;

  // take stellar parameters over to the binary
  bconfig.set(MB, BCO1) = ns1_config(MB);
  bconfig.set(MB, BCO2) = ns2_config(MB);

  bconfig.set(CHI, BCO1) = ns1_config(CHI);
  bconfig.set(CHI, BCO2) = ns2_config(CHI);

  bconfig.set(OMEGA, BCO1) = ns1_config(OMEGA);
  bconfig.set(OMEGA, BCO2) = ns2_config(OMEGA);

  bconfig.set(MADM, BCO1) = ns1_config(MADM);
  bconfig.set(MADM, BCO2) = ns2_config(MADM);

  // output configuration file
  bconfig.write_config(ss.str()+".info");

  // print final configuration
  std::cout << "Imported binary config:" << std::endl;
  std::cout << bconfig << std::endl;

  // create copy of data file, none of the fields are actually changed yet
  std::string data_fn = bin_fn.substr(0, bin_fn.size()-5)+".dat";
  fs::copy(data_fn, bconfig.config_outputdir()+ss.str()+".dat");
}
