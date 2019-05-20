#include <iostream>
#include <bitset>
#include <cmath>
#include <array>
#include <cassert>


using namespace std;


std::array<int, 1000*1000> checker{};

bool check(int i, int j, int t) {
    int& val = checker.at((500 + i) + 1000 * (500 + j));
    if (val == 0) {
        val = t;
        return true;
    } else {
        if (val + 1 == t) {
            val = t;
            return true;
        } else {
            return false;
        }
    }
}

volatile int unused = 0;

void print(int i, int j, int t) {
//  std::cout << i << " " << j << " " << t << endl;
//    cerr << "ijt: " << i << " " << j << " " << t << endl;
//  assert(check(i, j, t));
    unused += i + j + t;
}

int main() {
  std::array<int8_t,24> lut_ii{
      0, 1,-1, 0, 0, 1, 1,-1,-1, 1,-1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 1,-1
  };
  std::array<int8_t,24> lut_jj{
      0, 0, 0, 1,-1, 1,-1, 1,-1, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0
  };
  std::array<int8_t,24> lut_tt{
      0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1
  };
  std::array<int8_t,24> lut_mode {
      0,19,19,14,14, 0, 0, 0, 0,14,14,19,19, 0,14,14, 0,14,14,19,19, 0,19,19
  };

  int start_big = 0;

  int start_small_x = 14;

  int start_small_y = 19;

  // size of array is the logarithm of desired tile size
  std::array<int8_t, 10> state{};
  // the number of printed numbers higer than specified
  // by 10/7 = 1.42
  for (int x = 0; x < 70000000; x++) {
  //for (int x = 0; x < 4000; x++) {

    // step
    int32_t i = 0;
    int32_t j = 0;
    int32_t t = 0;

    for (size_t k = 0; k < state.size(); k++) {
      
      int32_t tt = lut_tt[state[k]];
      int32_t ii = lut_ii[state[k]];
      int32_t jj = lut_jj[state[k]];

      t += tt << k;
      i += ii << k;
      j += jj << k;

    }


    // print
    print(i, j, t);
    if (lut_mode[state[0]] == 0) {
      print(i, j, t + 1);
    }

    for (size_t k = 0; k < state.size(); k++) {
      if (state[k] != 13 && state[k] != 18 && state[k] != 23) {
        state[k]++;
        for (size_t kk = k; kk > 0; kk--) {
          int8_t mode = lut_mode[state[kk]];
          state[kk-1] = mode;
        }
        break;
      }
    }
  }
     
  cerr << "unused: " << unused << endl;
}
