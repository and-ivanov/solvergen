#ifndef SOLVERHELPER_H
#define SOLVERHELPER_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>

#ifndef NDEBUG
#include <fenv.h>
#endif

#include "vectorclass.h"
#include "immintrin.h"


using Int = int64_t;
using Real = double;

constexpr Real RealEps = std::numeric_limits<Real>::epsilon();

namespace SoverHelperImpl {

void check(const std::string& msg, bool expr) {
    if (!expr) {
        std::cerr << msg << std::endl;
        std::exit(1);
    }
}

}; // namespace

#define CHECK(expr) \
    SoverHelperImpl::check("Condition is not satisfied: " #expr, expr)

#define CHECK2(expr, msg) \
    SoverHelperImpl::check(msg, expr)

//***************Tiling methods************************//
enum TileKind : int8_t {
    TileLarge = 0, TileSmallX = 1, TileSmallY = 2
};

namespace TilesNamespace {

    constexpr const int8_t largeTileSize = 14;
    constexpr const int8_t smallTileSize = 5;
    constexpr const int8_t largeTileStart = 0;
    constexpr const int8_t smallTileXStart = largeTileStart + largeTileSize;
    constexpr const int8_t smallTileYStart = smallTileXStart + smallTileSize;
    constexpr const int8_t lutSize = smallTileYStart + smallTileSize;

    constexpr const std::array<int8_t, lutSize> lutI{0, 0, 0,-1, 1,-1, 1,-1, 1,-1, 1, 0, 0, 0,/**/-1, 1, 0, 0, 0,/**/ 0, 0, 0,-1, 1};
    constexpr const std::array<int8_t, lutSize> lutJ{0,-1, 1, 0, 0,-1,-1, 1, 1, 0, 0,-1, 1, 0,/**/ 0, 0, 0,-1, 1,/**/-1, 1, 0, 0, 0};
    constexpr const std::array<int8_t, lutSize> lutT{0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,/**/ 0, 0, 0, 1, 1,/**/ 0, 0, 0, 1, 1};
    constexpr const std::array<int8_t, lutSize> lutK{0, 1, 1, 2, 2, 0, 0, 0, 0, 1, 1, 2, 2, 0,/**/ 1, 1, 0, 1, 1,/**/ 2, 2, 0, 2, 2};
};

template <int Level, typename F>
struct Tile {
    __attribute__((always_inline))
    inline static void run(TileKind tk, int i, int j, int t, int nx, int ny, int nt, F func) {
        using namespace TilesNamespace;

        constexpr int tileSize = 1 << (Level - 1);
        constexpr int tileWidth = tileSize * 2 - 1;
        int tileHeight = (tk == TileLarge) ? (tileSize * 4 - 1) : (tileSize * 2 - 1);

        if (i + tileWidth < 0 || i - tileWidth >= nx ||
            j + tileWidth < 0 || j - tileWidth >= ny ||
            t + tileHeight < 0 || t >= nt) return;

        int8_t llstart = 0;
        int8_t llend = 0;

        switch (tk) {
            case TileLarge:
                llstart = largeTileStart;
                llend = llstart + largeTileSize;
                break;
            case TileSmallX:
                llstart = smallTileXStart;
                llend = llstart + smallTileSize;
                break;
            case TileSmallY:
                llstart = smallTileYStart;
                llend = llstart + smallTileSize;
                break;
        }

        for (int ll = llstart; ll < llend; ll++) {
            Tile<Level - 1, F>::run(
                    TileKind(lutK[ll]),
                    i + tileSize * lutI[ll],
                    j + tileSize * lutJ[ll],
                    t + tileSize * lutT[ll],
                    nx, ny, nt, func);
        }
    }
};

template <typename F> struct Tile<0, F> {
    __attribute__((always_inline))
    inline static void run(TileKind tk, int i, int j, int t, int nx, int ny, int nt, F func) {
        if (0 <= i && i < nx && 0 <= j && j < ny)
        {
            if (0 <= t && t < nt)
                func(i, j, t);
            if (tk == TileLarge)
                if (0 <= t + 1 && t + 1 < nt)
                    func(i, j, t + 1);
        }
    }
};

template <typename F>
void tiledLoopsTemplated(int nx, int ny, int nt, F f) {
    constexpr int tileSize = 15;
    int nodesWidth = std::max(nx, std::max(ny, nt));
    int requiredTileSize = 0;
    while ((1 << requiredTileSize) < nodesWidth) requiredTileSize++;
    assert(requiredTileSize <= tileSize);
    Tile<tileSize, decltype(f)>::run(TileLarge, 0, 0, -(1 << tileSize), nx, ny, nt, f);
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

    int base = std::max(nx, ny);
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

template <typename F>
void classicLoops(int nx, int ny, int nt, F func) {
    for (int t = 0; t < nt; t++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                func(i, j, t);
}


//***************End of tiling methods*****************//


struct Size {
    explicit Size(Int x, Int y) : x(x), y(y) {}
    Int x, y;
};

bool operator==(const Size& l, const Size& r) {
    return l.x == r.x && l.y == r.y;
}

struct Range {
    Range() {}

    Range& x(Int s, Int e) {
        xs = s;
        xe = e;
        return *this;
    }
    Range& y(Int s, Int e) {
        ys = s;
        ye = e;
        return *this;
    }
    Range dx(Int d) {
        return Range().x(xs + d, xe + d).y(ys, ye);
    }
    Range dy(Int d) {
        return Range().x(xs, xe).y(ys + d, ye + d);
    }

    Size size() const { return Size(xe - xs, ye - ys); }

    Int xs = 0;
    Int xe = -1;
    Int ys = 0;
    Int ye = -1;
};

// Modifialble: Temp | Var
// Leaf: Modifiable | Con
// Expr: Leaf | Leaf + Expr
// Assign: Modifiable = Expr
// Sequence: Assign | Assign , Sequence

template <typename T> struct isModifiable { static const bool value = false; };

template <typename T> struct isExprNode {
    static const bool value = std::is_arithmetic<T>::value;
};

template <typename T> struct isStrictExprNode {
    static const bool value = isExprNode<T>::value && !std::is_arithmetic<T>::value;
};

template <typename T> struct isSequence { static const bool value = false; };

template <typename T1, typename T2>
struct hasStrictExprNode {
    static const bool value = isExprNode<T1>::value && isExprNode<T1>::value &&
        (isStrictExprNode<T1>::value || isStrictExprNode<T2>::value);
};


struct Expr;
struct Var;

struct Assign {};

struct Sequence {};

struct Array {
    explicit Array(Size size) : mSize(size), mData(size_t(size.x * size.y)) {}
    explicit Array(Int nx, Int ny) : Array(Size(nx, ny)) {}

    Array& load(const std::string& s) {
        std::ifstream ifs(s);
        CHECK2(ifs.is_open(), "Can't open file for reading: " + s);
        ifs.seekg(0, std::ios_base::end);
        auto fileSize = ifs.tellg();
        ifs.seekg(0, std::ios_base::beg);
        Int arraySize = mSize.x * mSize.y;
        CHECK2(fileSize == arraySize * sizeof(Real),
               "Size of data in file '" + s + "' is not equal to requested" +
               "(" + std::to_string(mSize.x) +
               "," + std::to_string(mSize.y) + ")");
        ifs.read(reinterpret_cast<char*>(mData.data()), arraySize * sizeof(Real));
        return *this;
    }
    Array& save(const std::string& s) {
        std::ofstream ofs(s);
        CHECK2(ofs.is_open(), "Can't open file for writing: " + s);
        Int arraySize = mSize.x * mSize.y;
        ofs.write(reinterpret_cast<char*>(mData.data()), arraySize * sizeof(Real));
        return *this;
    }

    Real& val(Int x, Int y) __attribute__((always_inline,flatten)) {
        assert(0 <= x && x < mSize.x);
        assert(0 <= y && y < mSize.y);
        return mData[size_t(x + mSize.x * y)];
    }

    Assign operator=(const Array&) { return Assign(); }
    template <typename Expr, std::enable_if_t<isExprNode<Expr>::value, int> = 0>
    Assign operator=(const Expr&) { return Assign(); }

    const Size& size() const { return mSize; }
    const Range range() const { return Range().x(0, mSize.x).y(0, mSize.y); }

    Size mSize;
    std::vector<Real> mData;
};


struct Var {
    Var(Array& a, Range r = Range()) : mArray(a), mRange(r) {
        if (mRange.xe == -1) mRange.xe = a.size().x;
        if (mRange.ye == -1) mRange.ye = a.size().y;

        assert(0 <= mRange.xs);
        assert(mRange.xs <= mRange.xe);
        assert(mRange.xe <= mArray.size().x);

        assert(0 <= mRange.ys);
        assert(mRange.ys <= mRange.ye);
        assert(mRange.ye <= mArray.size().y);
    }
    Var(const Var& var) = default;

    Var& x(Int s, Int e) { mRange.xs = s; mRange.xe = e; return *this; }
    Var& y(Int s, Int e) { mRange.ys = s; mRange.ye = e; return *this; }

    Var dx(Int d) { return Var(mArray, mRange.dx(d)); }
    Var dy(Int d) { return Var(mArray, mRange.dy(d)); }

    Assign operator=(const Var&) { return Assign(); }
    template <typename Expr, std::enable_if_t<isExprNode<Expr>::value, int> = 0>
    Assign operator=(const Expr&) { return Assign(); }

    Assign operator+=(const Expr&) { return Assign(); }
    Assign operator-=(const Expr&) { return Assign(); }
    Assign operator*=(const Expr&) { return Assign(); }
    Assign operator/=(const Expr&) { return Assign(); }

    __attribute__((always_inline,flatten))
    inline Real& val(Int x, Int y) __attribute__((always_inline,flatten)) {
        assert(0 <= x && x < mRange.size().x);
        assert(0 <= y && y < mRange.size().y);
        return mArray.val(mRange.xs + x, mRange.ys + y);
    }

    const Range& range() const { return mRange; }

    Array& mArray;
    Range mRange;
};

struct Temp {
    Temp() {}

    Assign operator=(const Temp&) { return Assign(); }
    template <typename Expr, std::enable_if_t<isExprNode<Expr>::value, int> = 0>
    Assign operator=(const Expr&) { return Assign(); }

    Assign operator+=(const Expr&) { return Assign(); }
    Assign operator-=(const Expr&) { return Assign(); }
    Assign operator*=(const Expr&) { return Assign(); }
    Assign operator/=(const Expr&) { return Assign(); }
};

template <Int Dim>
struct Idx {};

template <Int Dim>
struct Nodes {};

template <> struct isModifiable<Temp> { static const bool value = true; };
template <> struct isModifiable<Var> { static const bool value = true; };
template <> struct isModifiable<Array> { static const bool value = true; };


template <> struct isExprNode<Array> { static const bool value = true; };
template <> struct isExprNode<Expr> { static const bool value = true; };
template <> struct isExprNode<Temp> { static const bool value = true; };
template <> struct isExprNode<Var> { static const bool value = true; };
template <Int Dim> struct isExprNode<Idx<Dim>> { static const bool value = true; };
template <Int Dim> struct isExprNode<Nodes<Dim>> { static const bool value = true; };

template <> struct isSequence<Assign> { static const bool value = true; };
template <> struct isSequence<Sequence> { static const bool value = true; };


struct Expr {
    Expr() {}
    template <typename T, std::enable_if_t<isExprNode<Expr>::value, int> = 0>
    Expr(const T&) {}
};

template <typename Expr1, typename Expr2,
          std::enable_if_t<hasStrictExprNode<Expr1, Expr2>::value, int> = 0>
Expr max(Expr1, Expr2) { return Expr(); }


template <typename Expr1, typename Expr2,
          std::enable_if_t<hasStrictExprNode<Expr1, Expr2>::value, int> = 0>
Expr min(Expr1, Expr2) { return Expr(); }

template <typename Expr1, typename Expr2,
          std::enable_if_t<hasStrictExprNode<Expr1, Expr2>::value, int> = 0>
Expr operator+(Expr1, Expr2) { return Expr(); }

template <typename Expr1, typename Expr2,
          std::enable_if_t<hasStrictExprNode<Expr1, Expr2>::value, int> = 0>
Expr operator-(Expr1, Expr2) { return Expr(); }

template <typename Expr1, typename Expr2,
          std::enable_if_t<hasStrictExprNode<Expr1, Expr2>::value, int> = 0>
Expr operator*(Expr1, Expr2) { return Expr(); }

template <typename Expr1, typename Expr2,
          std::enable_if_t<hasStrictExprNode<Expr1, Expr2>::value, int> = 0>
Expr operator/(Expr1, Expr2) { return Expr(); }

template <typename E, std::enable_if_t<isStrictExprNode<E>::value, int> = 0>
Expr operator+(E) { return Expr(); }

template <typename E, std::enable_if_t<isStrictExprNode<E>::value, int> = 0>
Expr operator-(E) { return Expr(); }

template <typename Expr1, typename Expr2,
          std::enable_if_t<hasStrictExprNode<Expr1, Expr2>::value, int> = 0>
Expr operator+=(Expr1, Expr2) { return Expr(); }


template <typename Expr1, typename Expr2,
          std::enable_if_t<hasStrictExprNode<Expr1, Expr2>::value, int> = 0>
bool operator>(Expr1, Expr2) { return true; }

template <typename Expr1, typename Expr2,
          std::enable_if_t<hasStrictExprNode<Expr1, Expr2>::value, int> = 0>
bool operator<(Expr1, Expr2) { return true; }

template <typename S1, typename S2,
          std::enable_if_t<isSequence<S1>::value, int> = 0,
          std::enable_if_t<isSequence<S2>::value, int> = 0>
Sequence operator,(S1, S2) { return Sequence(); }


Expr pow(Expr, Expr) { return Expr(); }
Expr exp(Expr) { return Expr(); }


inline std::string withLeadingZeros(Int value, Int leadingZeros)
{
     std::ostringstream oss;
     oss << std::setw(leadingZeros) << std::setfill('0') << value;
     return oss.str();
}

// opearators disambiguation

Real min(Real lhs, Real rhs) {
    return std::min(lhs, rhs);
}

Real max(Real lhs, Real rhs) {
    return std::max(lhs, rhs);
}

// vectorization

struct Real4 {
    Real4() {}
    Real4(Vec4d val) : val(val) {}

    Real4& operator+() {
        return *this;
    }
    const Real4& operator-() {
        val = -val;
        return *this;
    }

    const Real4& operator=(Real rhs) {
        val = rhs;
        return *this;
    }

    const Real4& operator+=(const Real4& rhs) {
        val += rhs.val;
        return *this;
    }
    const Real4& operator-=(const Real4& rhs) {
        val -= rhs.val;
        return *this;
    }
    const Real4& operator*=(const Real4& rhs) {
        val *= rhs.val;
        return *this;
    }
    const Real4& operator/=(const Real4& rhs) {
        val /= rhs.val;
        return *this;
    }

    const Real4& operator+=(Real rhs) {
        val += rhs;
        return *this;
    }
    const Real4& operator-=(Real rhs) {
        val -= rhs;
        return *this;
    }
    const Real4& operator*=(Real rhs) {
        val *= rhs;
        return *this;
    }
    const Real4& operator/=(Real rhs) {
        val /= rhs;
        return *this;
    }

    Vec4d val;
};

class DefineLoop {
public:
    using DefineLoopFn = std::function<void(int, int)>;
    using rhsFn = std::function<void()>;
    DefineLoop(DefineLoopFn)  {}
    DefineLoop() {}
    DefineLoop operator+(rhsFn) { return DefineLoop(); }

};

class LoopSequence {
public:
    LoopSequence() {}
};

//************************(lambaNum, vecN, i, j)***
LoopSequence initLoopSequence(int, int, int, int) {
    return LoopSequence();
}

Real4 operator+(const Real4& lhs, const Real4& rhs) {
    return Real4(lhs.val + rhs.val);
}

Real4 operator-(const Real4& lhs, const Real4& rhs) {
    return Real4(lhs.val - rhs.val);
}

Real4 operator*(const Real4& lhs, const Real4& rhs) {
    return Real4(lhs.val * rhs.val);
}

Real4 operator/(const Real4& lhs, const Real4& rhs) {
    return Real4(lhs.val / rhs.val);
}

Real4 operator+(const Real& lhs, const Real4& rhs) {
    return Real4(lhs + rhs.val);
}

Real4 operator+(const Real4& lhs, const Real& rhs) {
    return Real4(lhs.val + rhs);
}

Real4 operator-(const Real& lhs, const Real4& rhs) {
    return Real4(lhs - rhs.val);
}

Real4 operator-(const Real4& lhs, const Real& rhs) {
    return Real4(lhs.val - rhs);
}

Real4 operator*(const Real& lhs, const Real4& rhs) {
    return Real4(lhs * rhs.val);
}

Real4 operator*(const Real4& lhs, const Real& rhs) {
    return Real4(lhs.val * rhs);
}

Real4 operator/(const Real& lhs, const Real4& rhs) {
    return Real4(lhs / rhs.val);
}

Real4 operator/(const Real4& lhs, const Real& rhs) {
    return Real4(lhs.val / rhs);
}

Real4 min(const Real& lhs, const Real4& rhs) {
    Vec4d l(lhs);
    return Real4(min(l, rhs.val));
}

Real4 min(const Real4& lhs, const Real& rhs) {
    Vec4d r(rhs);
    return Real4(min(lhs.val, r));
}

Real4 min(const Real4& lhs, const Real4& rhs) {
    return Real4(min(lhs.val, rhs.val));
}

Real4 max(const Real& lhs, const Real4& rhs) {
    Vec4d l(lhs);
    return Real4(max(l, rhs.val));
}

Real4 max(const Real4& lhs, const Real& rhs) {
    Vec4d r(rhs);
    return Real4(max(lhs.val, r));
}

Real4 max(const Real4& lhs, const Real4& rhs) {
    return Real4(max(lhs.val, rhs.val));
}

__attribute__((always_inline))
inline void loadFromPtr(Real& l, Real* r) {
    l = *r;
}

__attribute__((always_inline))
inline void loadFromPtr(Real4& l, Real* r) {
    l.val.load(r);
}

__attribute__((always_inline))
inline void storeToPtr(Real* l, Real r) {
    *l = r;
}

__attribute__((always_inline))
inline void storeToPtr(Real* l, Real4 r) {
    r.val.store(l);
}

#endif // SOLVERHELPER_H
