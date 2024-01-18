#include <vector>
#include <iostream>
template<class reader_t>
void interp_data(reader_t& input_reader, std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& zz) {
  using namespace std;
  auto const Npts = xx.size();
  std::vector<typename reader_t::output_ary_t> all_data(Npts);
  #pragma omp parallel for firstprivate(input_reader)
  for(auto i = 0; i < Npts; ++i) {
    all_data[i] = input_reader.export_pointwise(xx[i], yy[i], zz[i]);
  }
  
  for(auto i = 0; i < Npts; ++i) {
    std::cout << "(" << xx[i] << ", " << yy[i] << ", " << zz[i] << ") - ";
    
    for(auto& e : all_data[i])
      cout << e << " - ";
    cout << endl;
  }
  cout << xx[200] << ", " << yy[200] << ", " << zz[200] << ", " 
    << all_data[200][reader_t::OUTPUT_VARS::ALPHA] << "\n\t"
    << all_data[200][reader_t::OUTPUT_VARS::K12] << ", "
    << all_data[200][reader_t::OUTPUT_VARS::K13] << ", "
    << all_data[200][reader_t::OUTPUT_VARS::K23] << "\n";

}