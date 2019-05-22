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

template <typename F>
void largeTile(int ox, int oy, int ot,
        int nx, int ny, int nt, int smallTileHeight, F func) 
{
    int largeTileHeight = smallTileHeight * 2;
    for (int t = max(0, ot); t < min(nt, ot + largeTileHeight); t++) {
        int at = t - ot;
        int halfWidth = std::min(at, largeTileHeight - at - 1);
        for (int j = max(0, oy - halfWidth); j <= min(ny - 1, oy + halfWidth); j++) {
            for (int i = max(0, ox - halfWidth); i <= min(nx - 1, ox + halfWidth); i++) {
                func(i, j, t);
            }
        }
    }
}

template <typename F>
void smallTileX(int ox, int oy, int ot,
        int nx, int ny, int nt, int smallTileHeight, F func) 
{
    int largeTileHeight = smallTileHeight * 2;
    //cerr << "small tile x" << endl;
    for (int t = max(0, ot); t < min(nt, ot + smallTileHeight); t++) {
        //cerr << "t: " << t << endl;
        int at = t - ot;
        //cerr << "at: " << at << endl;
        int halfWidthX = smallTileHeight - at - 1;
        int halfWidthY = at;
        //cerr << "halfWidthX: " << halfWidthX << endl;
        //cerr << "halfWidthY: " << halfWidthY << endl;
        for (int j = max(0, oy - halfWidthY); j <= min(ny - 1, oy + halfWidthY); j++) {
            for (int i = max(0, ox - halfWidthX); i <= min(nx - 1, ox + halfWidthX); i++) {
                func(i, j, t);
            }
        }
    }
}

template <typename F>
void smallTileY(int ox, int oy, int ot,
        int nx, int ny, int nt, int smallTileHeight, F func) 
{
    int largeTileHeight = smallTileHeight * 2;
    for (int t = max(0, ot); t < min(nt, ot + smallTileHeight); t++) {
        int at = t - ot;
        int halfWidthX = at;
        int halfWidthY = smallTileHeight - at - 1;
        for (int j = max(0, oy - halfWidthY); j <= min(ny - 1, oy + halfWidthY); j++) {
            for (int i = max(0, ox - halfWidthX); i <= min(nx - 1, ox + halfWidthX); i++) {
                func(i, j, t);
            }
        }
    }
}

template <typename F>
void fixedTiledLoops(int nx, int ny, int nt, int smallTileHeight, F func) {
    int largeTileHeight = smallTileHeight * 2; 
    int ntt = (nt - 1) / largeTileHeight + 1;
    int nxx = (nx - 1) / largeTileHeight + 1;
    int nyy = (ny - 1) / largeTileHeight + 1;

    // because tiles will not cover the last layers
    ntt++;
    nxx++;
    nyy++;

    for (int tt = 0; tt < ntt; tt++) {
        int ot = tt * largeTileHeight;

        //cerr << "L1" << endl;
        for (int yy = 0; yy < nyy; yy++) {
            int oy = yy * largeTileHeight;
            for (int xx = 0; xx < nxx; xx++) {
                int ox = xx * largeTileHeight;
                largeTile(
                        ox + smallTileHeight,
                        oy + smallTileHeight,
                        ot - smallTileHeight, 
                        nx, ny, nt, smallTileHeight, func);
            }
        }
        //cerr << "X1" << endl;
        for (int yy = 0; yy < nyy; yy++) {
            int oy = yy * largeTileHeight;
            for (int xx = 0; xx < nxx; xx++) {
                int ox = xx * largeTileHeight;
                smallTileX(
                        ox + smallTileHeight,
                        oy,
                        ot, 
                        nx, ny, nt, smallTileHeight, func);
            }
        }
        //cerr << "Y1" << endl;
        for (int yy = 0; yy < nyy; yy++) {
            int oy = yy * largeTileHeight;
            for (int xx = 0; xx < nxx; xx++) {
                int ox = xx * largeTileHeight;
                smallTileY(
                        ox,
                        oy + smallTileHeight,
                        ot, 
                        nx, ny, nt, smallTileHeight, func);
            }
        }
        //cerr << "L2" << endl;
        for (int yy = 0; yy < nyy; yy++) {
            int oy = yy * largeTileHeight;
            for (int xx = 0; xx < nxx; xx++) {
                int ox = xx * largeTileHeight;
                largeTile(
                        ox,
                        oy,
                        ot, 
                        nx, ny, nt, smallTileHeight, func);
            }
        }
        //cerr << "X2" << endl;
        for (int yy = 0; yy < nyy; yy++) {
            int oy = yy * largeTileHeight;
            for (int xx = 0; xx < nxx; xx++) {
                int ox = xx * largeTileHeight;
                smallTileX(
                        ox,
                        oy + smallTileHeight,
                        ot + smallTileHeight, 
                        nx, ny, nt, smallTileHeight, func);
            }
        }
        //cerr << "X2" << endl;
        for (int yy = 0; yy < nyy; yy++) {
            int oy = yy * largeTileHeight;
            for (int xx = 0; xx < nxx; xx++) {
                int ox = xx * largeTileHeight;
                smallTileY(
                        ox + smallTileHeight,
                        oy,
                        ot + smallTileHeight, 
                        nx, ny, nt, smallTileHeight, func);
            }
        }
    }


}

int main() {
    int nx = 400;
    int ny = 400;
    int nt = 100;
    
    volatile int unused = 0;

    auto f = [&unused](int i, int j, int t) {
//        std::cout << i << " " << j << " " << t << endl;
        //cerr << "ijt: " << i << " " << j << " " << t << endl;
        assert(check(i, j, t));
        unused += i + j + t;
    };
    fixedTiledLoops(nx, ny, nt, 5, f);

    cerr << "unused: " << unused << endl;
}
