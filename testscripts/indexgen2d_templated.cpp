#include <iostream>
#include <vector>
#include <bitset>
#include <cmath>
#include <array>
#include <cassert>


using namespace std;


class Checker {
public:
    Checker(int nx, int ny)
        : nx(nx), ny(ny), data(nx * ny, -1)
    {
    }
    bool check(int i, int j, int t) {
        int& val = data.at(i + nx * j);
        auto is = {-1, 0, 1};
        auto js = {-1, 0, 1};
        for (int iss : is) {
            for (int jss : js) {
                int ii = i + iss;
                int jj = j + jss;
                if (ii < 0 || jj < 0 || ii >= nx || jj >= ny) continue;
                int vval = data.at(ii + nx * jj);
                if (vval < t - 1) {
                    return false;
                }
            }
        }
        if (val == t - 1) {
            val = t;
            return true;
        } else {
            return false;
        }
    }
private:    
    int nx, ny;
    std::vector<int> data;
};


template <typename T>
inline T pow2(T k) {
    return 1ll << k;
}

template <typename F>
void tiledLoops(int nx, int ny, int nt, F func) {
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
    
    int base = max(nx, ny);
    base = base | 1; // round up to the nearest odd number
    int height = nt;
    
    int32_t tileParam = 1;
    while (pow2(tileParam+1) < base + height - 1) {
        tileParam++;
    }
    assert(pow2(tileParam+1) >= base + height - 1);
    //cerr << "tileParam: " << tileParam << endl;
    
    int32_t sx = - base / 2;
    //cerr << "sx: " << sx << endl;
    int32_t ex = sx + nx - 1;
    //cerr << "ex: " << ex << endl;
    
    int32_t sy = - base / 2;
    //cerr << "sy: " << sy << endl;
    int32_t ey = sy + ny - 1;
    //cerr << "ey: " << ey << endl;
    
    int32_t st = + base / 2;
    //cerr << "st: " << st << endl;
    int32_t et = st + nt - 1;
    //cerr << "et: " << et << endl;
    
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
    //cerr << "iterations: " << iterations << endl;
    
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
                    func(iis[0] - sx, jjs[0] - sy, tts[0] - st);
                }
                if (lut_mode[state[0]] == 0) {
                    if(st <= tts[0] + 1 && tts[0] + 1 <= et) {
                        func(iis[0] - sx, jjs[0] - sy, tts[0] + 1 - st);
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
     
}

template <int Level, typename F>
struct SmallTileX {
    static void run(int i, int j, int t, int nx, int ny, int nt, F func);
};
template <int Level, typename F>
struct SmallTileY {
    static void run(int i, int j, int t, int nx, int ny, int nt, F func);
};
template <int Level, typename F>
struct LargeTile {
    static void run(int i, int j, int t, int nx, int ny, int nt, F func);
};

template <int Level, typename F>
void SmallTileX<Level, F>::run(int i, int j, int t, int nx, int ny, int nt, F func) {
    constexpr int tileSize = 1 << (Level - 1);
    constexpr int tileWidth = tileSize * 2 - 1;
    constexpr int tileHeight = tileSize * 2 - 1;

    if (i + tileWidth < 0 || i - tileWidth >= nx ||
        j + tileWidth < 0 || j - tileWidth >= ny ||
        t + tileHeight < 0 || t >= nt) return;

    SmallTileX<Level - 1, F>::run(i - tileSize, j, t, nx, ny, nt, func);
    SmallTileX<Level - 1, F>::run(i + tileSize, j, t, nx, ny, nt, func);

    LargeTile<Level - 1, F>::run(i, j, t, nx, ny, nt, func);

    SmallTileX<Level - 1, F>::run(i, j - tileSize, t + tileSize, nx, ny, nt, func);
    SmallTileX<Level - 1, F>::run(i, j + tileSize, t + tileSize, nx, ny, nt, func);
}

template <int Level, typename F>
void SmallTileY<Level, F>::run(int i, int j, int t, int nx, int ny, int nt, F func) {
    constexpr int tileSize = 1 << (Level - 1);
    constexpr int tileWidth = tileSize * 2 - 1;
    constexpr int tileHeight = tileSize * 2 - 1;

    if (i + tileWidth < 0 || i - tileWidth >= nx ||
        j + tileWidth < 0 || j - tileWidth >= ny ||
        t + tileHeight < 0 || t >= nt) return;

    SmallTileY<Level - 1, F>::run(i, j - tileSize, t, nx, ny, nt, func);
    SmallTileY<Level - 1, F>::run(i, j + tileSize, t, nx, ny, nt, func);

    LargeTile<Level - 1, F>::run(i, j, t, nx, ny, nt, func);

    SmallTileY<Level - 1, F>::run(i - tileSize, j, t + tileSize, nx, ny, nt, func);
    SmallTileY<Level - 1, F>::run(i + tileSize, j, t + tileSize, nx, ny, nt, func);
}

template <int Level, typename F>
void LargeTile<Level, F>::run(int i, int j, int t, int nx, int ny, int nt, F func) {
    constexpr int tileSize = 1 << (Level - 1);
    constexpr int tileWidth = tileSize * 2 - 1;
    constexpr int tileHeight = tileSize * 4 - 1;

    if (i + tileWidth < 0 || i - tileWidth >= nx ||
        j + tileWidth < 0 || j - tileWidth >= ny ||
        t + tileHeight < 0 || t >= nt) return;

    LargeTile<Level - 1, F>::run(i, j, t, nx, ny, nt, func);

    SmallTileX<Level - 1, F>::run(i, j - tileSize, t + tileSize, nx, ny, nt, func);
    SmallTileX<Level - 1, F>::run(i, j + tileSize, t + tileSize, nx, ny, nt, func);
    SmallTileY<Level - 1, F>::run(i - tileSize, j, t + tileSize, nx, ny, nt, func);
    SmallTileY<Level - 1, F>::run(i + tileSize, j, t + tileSize, nx, ny, nt, func);

    LargeTile<Level - 1, F>::run(i - tileSize, j - tileSize, t + tileSize, nx, ny, nt, func);
    LargeTile<Level - 1, F>::run(i + tileSize, j - tileSize, t + tileSize, nx, ny, nt, func);
    LargeTile<Level - 1, F>::run(i - tileSize, j + tileSize, t + tileSize, nx, ny, nt, func);
    LargeTile<Level - 1, F>::run(i + tileSize, j + tileSize, t + tileSize, nx, ny, nt, func);

    SmallTileX<Level - 1, F>::run(i - tileSize, j, t + 2 * tileSize, nx, ny, nt, func);
    SmallTileX<Level - 1, F>::run(i + tileSize, j, t + 2 * tileSize, nx, ny, nt, func);
    SmallTileY<Level - 1, F>::run(i, j - tileSize, t + 2 * tileSize, nx, ny, nt, func);
    SmallTileY<Level - 1, F>::run(i, j + tileSize, t + 2 * tileSize, nx, ny, nt, func);
    
    LargeTile<Level - 1, F>::run(i, j, t + 2 * tileSize, nx, ny, nt, func);
}

template <typename F> struct SmallTileX<0, F> {
    static void run(int i, int j, int t, int nx, int ny, int nt, F func) {
        if (0 <= i && i < nx &&
            0 <= j && j < ny &&
            0 <= t && t < nt) 
        {
            func(i, j, t);
        }
    }
};
template <typename F> struct SmallTileY<0, F> {
    static void run(int i, int j, int t, int nx, int ny, int nt, F func) {
        if (0 <= i && i < nx &&
            0 <= j && j < ny &&
            0 <= t && t < nt) 
        {
            func(i, j, t);
        }
    }
};
template <typename F> struct LargeTile<0, F> {
    static void run(int i, int j, int t, int nx, int ny, int nt, F func) {
        if (0 <= i && i < nx &&
            0 <= j && j < ny) 
        {
            if (0 <= t && t < nt) func(i, j, t);
            if (0 <= t + 1 && t + 1 < nt) func(i, j, t + 1);
        }
    }
};

int main() {
    int nx = 2000;
    int ny = 2000;
    int nt = 100;
    
    Checker checker(nx, ny);

    int unused = 0;

    //auto f = [&unused,&checker](int i, int j, int t) {
    //    std::cout << i << " " << j << " " << t << endl;
    //    //cerr << "ijt: " << i << " " << j << " " << t << endl;
    //    assert(checker.check(i, j, t));
    //    unused += i + j + t;
    //};
    //tiledLoops(nx, ny, nt, f);
    //templatedTiledLoops();
    auto f = [&](int i, int j, int t) {
        //std::cerr << i << " " << j << " " << t << endl;
        //std::cout << i << " " << j << " " << t << endl;
        //assert(checker.check(i, j, t));
        unused += i + j + t;
    };
    constexpr int tileSize = 17;
    LargeTile<tileSize, decltype(f)>::run(0, 0, -1 << tileSize, nx, ny, nt, f);

    cerr << "unused: " << unused << endl;
}
