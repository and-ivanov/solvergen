#include <iostream>
#include <vector>
#include <bitset>
#include <cmath>
#include <array>
#include <cassert>


using namespace std;


std::array<int, 1000*1000> checker{};

bool check(int i, int j, int t) {
    int& val = checker.at((500 + i) + 1000 * (500 + j));
    array<int, 3> is = {-1, 0, 1};
    array<int, 3> js = {-1, 0, 1};
    if (val == 0) {
        val = t;
        return true;
    } else {
        for (int iss : is) {
            for (int jss : js) {
                int vval = checker.at((500 + i + iss) + 1000 * (500 + j + jss));
                if (vval != 0) {
                    if (vval < t - 1) return false;
                }
            }
        }
        if (val == 0 || val == t - 1) {
            val = t;
            return true;
        } else {
            return false;
        }
    }
}

volatile int unused = 0;

void print(int i, int j, int t) {
    std::cout << i << " " << j << " " << t << endl;
//    cerr << "ijt: " << i << " " << j << " " << t << endl;
    assert(check(i, j, t));
    unused += i + j + t;
}

template <typename T>
inline T pow2(T k) {
    return 1ll << k;
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
    
    int nx = 50;
    int ny = 40;
    int nt = 30;
    
    int base = max(nx, ny);
    base = base | 1; // round up to the nearest odd number
    int height = nt;
    
    int32_t tileParam = 1;
    while (pow2(tileParam+1) < base + height - 1) {
        tileParam++;
    }
    assert(pow2(tileParam+1) >= base + height - 1);
    cerr << "tileParam: " << tileParam << endl;
    
    int32_t sx = - base / 2;
    cerr << "sx: " << sx << endl;
    int32_t ex = sx + nx - 1;
    cerr << "ex: " << ex << endl;
    
    int32_t sy = - base / 2;
    cerr << "sy: " << sy << endl;
    int32_t ey = sy + ny - 1;
    cerr << "ey: " << ey << endl;
    
    int32_t st = + base / 2;
    cerr << "st: " << st << endl;
    int32_t et = st + nt - 1;
    cerr << "et: " << et << endl;
    
    // size of array is the logarithm of desired tile size
    std::vector<int8_t> state(tileParam);
    
    // height of the tile depends on the number of iterations:
    // num of iterations: 2^{3i+1}-2^i
    // i from 1 to inf
    // height big: 2^{i+1} ( range: [0;2^{i+1}-1] )
    // height small: 2^i ( range: [0;2^i-1] )
    // width: 2^{i+1} - 1 ( range: [-(2^i - 1); +(2^i - 1)] )
    int64_t iterations = pow2(3ll * tileParam + 1ll) - pow2(tileParam);
    //iterations = 20;
    cerr << "iterations: " << iterations << endl;
    
    std::vector<int32_t> tts(tileParam + 1);
    std::vector<int32_t> iis(tileParam + 1);
    std::vector<int32_t> jjs(tileParam + 1);
    
    size_t K = state.size() - 1;
    
    bool finished = false;
    while (1) {

        //cerr << "=====" << endl;
        // step

        //for (int i = 0; i < state.size(); i++) {
        //    cerr << int(state[i]) << " ";
        //}
        //cerr << endl;

        bool skipTile = false;
        while (1) {
            int32_t ss = state[K];

            int32_t tt = lut_tt[ss];
            int32_t ii = lut_ii[ss];
            int32_t jj = lut_jj[ss];

            tts[K] = tts[K+1] + (tt << K);
            iis[K] = iis[K+1] + (ii << K);
            jjs[K] = jjs[K+1] + (jj << K);

            int32_t min_t = tts[K] + 0;
            int32_t min_x = iis[K] - (pow2(K + 1) - 1);
            int32_t min_y = jjs[K] - (pow2(K + 1) - 1);
            
            int32_t mode = lut_mode[ss];
            int32_t height = pow2(K + (mode == 0 ? 2 : 1)) - 1;

            int32_t max_t = tts[K] + height - 1;
            int32_t max_x = iis[K] + (pow2(K + 1) - 1);
            int32_t max_y = jjs[K] + (pow2(K + 1) - 1);

            //cerr 
            //    << min_x << " "
            //    << max_x << " "
            //    << min_y << " "
            //    << max_y << " "
            //    << min_t << " "
            //    << max_t << " "
            //    << endl;
            
            if (max_t < st || min_t > et ||
                max_x < sx || min_x > ex || 
                max_y < sy || min_y > ey) 
            {
                skipTile = true;
                break;
            }

            if (K == 0) break; 

            state[K-1] = lut_mode[state[K]];
            
            K--;
        };

        //cerr << "skipTile: " << skipTile << endl;
        //cerr << "ijt: " << iis[0] << " " << jjs[0] << " " << tts[0] << endl;

        if (!skipTile) {

            // print
            if (sx <= iis[0] && iis[0] <= ex &&
                sy <= jjs[0] && jjs[0] <= ey) {
                if (st <= tts[0] && tts[0] <= et) {
                    print(iis[0] - sx, jjs[0] - sy, tts[0] - st);
                }
                if (lut_mode[state[0]] == 0) {
                    if(st <= tts[0] + 1 && tts[0] + 1 <= et) {
                        print(iis[0] - sx, jjs[0] - sy, tts[0] + 1 - st);
                    }
                }
            }

        } 

        while (state[K] == 13 || state[K] == 18 || state[K] == 23) {
            K++;
            if (K == state.size()) { finished = true; break; }
        }
        if (finished == true) break;

        state[K]++;


    }
     
    cerr << "unused: " << unused << endl;
}
