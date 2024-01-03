/*
    Copyright 2019 Ludwig Jens Papenfort

    This file is part of Kadath.

    Kadath is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Kadath is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kadath.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MEMORY_HH
#define MEMORY_HH

#include <iostream>
#include <cstdlib>
#include <unordered_map>
#include <map>
#include <queue>
#include <vector>
#include <memory>
#include <algorithm>

// only do this if really necessary, e.g. using an Intel compiler not capable of compiling the flat_hash_map below
// this is slower than the hash map
//#define KADATH_VECTORMAP
#define DEFAULT_KAD_MEM
#ifndef KADATH_VECTORMAP
#include "implementation/flat_hash_map.hpp"
#endif

namespace Kadath {
/*
  These classes define a thin wrapper around memory allocation calls,
  which are abundant throughout the Kadath library.
  By redirecting these calls and keeping the allocated memory by
  circumventing deallocations, significant runtime can be saved.

  This uses a fast hash map internally, to retrieve and store
  allocated memory chunks.

  In case a compiler is not able to compile the hash map,
  a fallback solution by a search over a vector is implemented.
*/


class coef_mem {
  private:
    static std::array<std::unique_ptr<double>,2> mem_ptrs;
    static std::array<size_t,2> lengths;

  public:
    static double* get_mem(size_t const num, size_t const len) {
      #ifndef DEFAULT_KAD_MEM
      if(len > lengths[num]) {
      	mem_ptrs[num].reset(new double[len]);

        lengths[num] = len;
      }
      return mem_ptrs[num].get();
      #else
      return new double [len];
      #endif
    }
};

class MemoryMapper {
  private:
    struct FreeDeleter
    {
        void operator()(void* raw_mem_ptr) { std::free(raw_mem_ptr); }
    };

    using ptr_vec_t = std::vector<void*>;

    #ifdef KADATH_VECTORMAP

    using mem_map_t = std::vector<ptr_vec_t>;
    using mem_sizes_t = std::vector<size_t>;
    static mem_sizes_t mem_sizes;

    #else

    using mem_map_t = ska::flat_hash_map<size_t, ptr_vec_t>;

    #endif

    using ptr_t = std::unique_ptr<void,FreeDeleter>;
    using ptr_list_t = std::vector<ptr_t>;

    static mem_map_t mem_map;
    static ptr_list_t ptr_list;

  public:

    #ifndef DEFAULT_KAD_MEM
    static void* get_memory(size_t const sz) {

      if(sz == 0)
        return nullptr;

      void* raw_mem_ptr;

      #ifdef KADATH_VECTORMAP

      auto pos = std::find(std::begin(mem_sizes), std::end(mem_sizes), sz);
      if(pos == mem_sizes.end()) {
        mem_sizes.push_back(sz);
        mem_map.resize(mem_sizes.size());
        pos = --mem_sizes.end();
      }

      size_t index = pos - mem_sizes.begin();
      auto& mem = mem_map[index];

      #else

      auto& mem = mem_map[sz];

      #endif

      if(mem.empty()) {
        ptr_list.emplace_back(std::malloc(sz));
        raw_mem_ptr = ptr_list.back().get();
      }
      else {
        raw_mem_ptr = mem.back();
        mem.pop_back();
      }

      return raw_mem_ptr;
    }

    template<typename T>
    static T* get_memory(size_t const sz) {
      return static_cast<T*>(get_memory(sz * sizeof(T)));
    }

    static void release_memory(void* raw_mem_ptr, size_t const sz)  {
      if(raw_mem_ptr == nullptr)
        return;

      #ifdef KADATH_VECTORMAP
      auto pos = std::find(std::begin(mem_sizes), std::end(mem_sizes), sz);
      size_t index = pos - mem_sizes.begin();

      mem_map[index].push_back(raw_mem_ptr);
      #else
      mem_map[sz].push_back(raw_mem_ptr);
      #endif
    }

    template<typename T>
    static void release_memory(void* raw_mem_ptr, size_t const sz) {
      release_memory(raw_mem_ptr, sz * sizeof(T));
    }
    #else

    template<typename T>
    static T* get_memory(size_t const sz) {
      return new T [sz];
    }
   
    template<typename T>
    static void release_memory(T* mem_ptr, size_t const sz) {
      delete mem_ptr ;
    }
    template<typename T, size_t ary_sz>
    static void release_memory(T (*mem_ptr)[ary_sz], size_t const sz) {
      for(auto i = 0; i < ary_sz; ++i) {
        delete mem_ptr[i];
      }
      delete [] mem_ptr ;
    }
    #endif
};

struct MemoryMappable {
  #ifndef DEFAULT_KAD_MEM
  void* operator new(size_t sz)
  {
    return MemoryMapper::get_memory(sz);
  }

  void operator delete(void* mem_ptr, size_t const sz)
  {
    MemoryMapper::release_memory(mem_ptr, sz);
  }

  void* operator new[](size_t sz)
  {
    return MemoryMapper::get_memory(sz);
  }

  void operator delete[](void* mem_ptr, size_t const sz)
  {
    MemoryMapper::release_memory(mem_ptr, sz);
  }
  #endif
};
}
#endif
