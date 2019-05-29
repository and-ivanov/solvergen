#include "solverhelper.h"
#include "hrdata/hrdata.h"

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


//classicLoops - 26.22(Vladislav)
//tiledLoopsTemplated - 25.57(Vladislav)

//***************End of tiling methods*****************//

void extend_ghost(Var a, Int gh, Int gnx, Int gny) {
    for (Int i = 0; i < gh; i++) {
        a.x(i,i+1) = + a.x(gh,gh+1);
        a.x(gnx-i-1,gnx-i) = + a.x(gnx-gh-1,gnx-gh);
    }
    for (Int j = 0; j < gh; j++) {
        a.y(j,j+1) = + a.y(gh,gh+1);
        a.y(gny-j-1,gny-j) = + a.y(gny-gh-1,gny-gh);
    }
}

inline bool isNull(Expr a) { return -RealEps < a && a < RealEps; }

inline Expr limiter(Expr r) {
    return max(0.0, max(min(1.0, 2.0 * r), min(2.0, r)));
}

inline Expr advection(Expr c, Expr ppu, Expr pu, Expr u, Expr nu, Expr nnu) {
//    return - c * (u - pu);

    auto r1 = u - pu;
    auto r2 = nu - u;
    r1 += RealEps;
    r2 += RealEps;
    auto r = r1 / r2;

    auto pr1 = pu - ppu;
    auto pr2 = u - pu;
    pr1 += RealEps;
    pr2 += RealEps;
    auto pr = pr1 / pr2;

    auto nf12 = u + limiter(r) / 2.0 * (1.0 - c) * (nu - u);
    auto pf12 = pu + limiter(pr) / 2.0 * (1.0 - c) * (u - pu);

    return - c * (nf12 - pf12);
}


int main(int argc, char** argv) {
#ifndef NDEBUG
    feenableexcept(FE_INVALID | FE_OVERFLOW); // gnu extension
#endif

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " config.hrdata" << std::endl;
        std::exit(-1);
    }

    hrdata::Node conf = hrdata::parseFile(argv[1]);

    auto nx = conf["nx"].as<Int>(); CHECK(nx > 0);
    auto ny = conf["ny"].as<Int>(); CHECK(ny > 0);

    auto time_steps = conf["time_steps"].as<Int>(); CHECK(time_steps > 0);

    auto dx = conf["dx"].as<Real>(); CHECK(dx > 0);
    auto dy = conf["dy"].as<Real>(); CHECK(dy > 0);
    auto dt = conf["time_step"].as<Real>(); CHECK(dt > 0);

    Int gh = 2;

    Int gnx = nx + 2 * gh;
    Int gny = ny + 2 * gh;

    Int lx_gh = gh;
    Int ly_gh = gh;

    Int rx_gh = gnx - gh;
    Int ry_gh = gny - gh;

    auto save_rate = conf["save_rate"].as<Int>();
    auto save_path = conf["save_path"].as<std::string>();

    auto inside = Range().x(lx_gh, rx_gh).y(ly_gh, ry_gh);

    auto size = Size(gnx, gny);

    // last letter a denotes entire Array
    // last letter i of Var denotes inner area of Array without ghost

    Array vxa(size), vya(size), sxxa(size), syya(size), sxya(size);


    Array cpa(size), csa(size), rhoa(size), laa(size), mua(size), lma(size);
    Array icpa(size), icsa(size), ilaa(size), imua(size), ilma(size);

    Array w1a(size), w2a(size), w3a(size), w4a(size);

    Var w1i(w1a, inside);
    Var w2i(w2a, inside);
    Var w3i(w3a, inside);
    Var w4i(w4a, inside);

    Var vxi(vxa, inside);
    Var vyi(vya, inside);
    Var sxxi(sxxa, inside);
    Var syyi(syya, inside);
    Var sxyi(sxya, inside);

    Var cpi(cpa, inside);
    Var csi(csa, inside);
    Var rhoi(rhoa, inside);
    Var lai(laa, inside);
    Var mui(mua, inside);
    Var lmi(lma, inside);

    Var icpi(icpa, inside);
    Var icsi(icsa, inside);
    Var imui(imua, inside);
    Var ilmi(ilma, inside);

    {
        auto pi = Array(nx, ny).load(conf["initial_pressure"].as<std::string>());
        sxxi = pi;
        syyi = pi;
    }

    {
        auto cpf = Array(nx, ny).load(conf["cp"].as<std::string>());
        cpi = cpf;
    }

    {
        auto csf = Array(nx, ny).load(conf["cs"].as<std::string>());
        csi = csf;
    }

    {
        auto rhof = Array(nx, ny).load(conf["rho"].as<std::string>());
        rhoi = rhof;
    }

    // extend material values on the boundaries
    extend_ghost(cpa, gh, gnx, gny);
    extend_ghost(csa, gh, gnx, gny);
    extend_ghost(rhoa, gh, gnx, gny);

    {
        mua = rhoa * csa * csa,
        laa = rhoa * cpa * cpa - 2 * mua;

        lma = laa + 2. * mua;

        icpa = 1. / (2. * cpa);
        icsa = 1. / (2. * csa);
        imua = 1. / (2. * mua);
        ilma = 1. / (2. * lma);
    }

    for (Int t = 0; t < time_steps; t++) {

        Temp c1, c2;
        Temp dw1, dw2, dw3, dw4;

        DefineLoop([&](int nx, int ny){
            Int xx = (nx / 4);
            Int yy = (ny / 2);
            Int tt = (time_steps * 4);
            tiledLoopsTemplated(xx, yy, tt,
                       [&](int i, int j, int t) __attribute__((always_inline)) {
                switch(t % 4) {
                case 0:
                {
                    initLoopSequence(0, 4, 4*i, 2*j);
                    initLoopSequence(0, 4, 4*i, 2*j+1);
                    break;
                }
                case 1:
                {
                    initLoopSequence(1, 4, 4*i, 2*j+1);
                    initLoopSequence(1, 4, 4*i, 2*j+1);
                    break;
                }
                case 2:
                {
                    initLoopSequence(2, 4, 4*i, 2*j+1);
                    initLoopSequence(2, 4, 4*i, 2*j+1);
                    break;
                }
                case 3:
                {
                    initLoopSequence(3, 4, 4*i, 2*j+1);
                    initLoopSequence(3, 4, 4*i, 2*j+1);
                    break;
                }
                }
            });
        })+
        [&]() { // to_omega_x
                    w1i = + vxi * icpi + sxxi * ilmi,
                    w2i = - vxi * icpi + sxxi * ilmi,
                    w3i = + vyi * icsi + sxyi * imui,
                    w4i = - vyi * icsi + sxyi * imui;
        }+
        [&]() { // omega_x_solve
                    c1 = cpi * dt / dx,
                    c2 = csi * dt / dx,

                    //        dw1 = advection(c1, w1i.dx(+2), w1i.dx(+1), w1i, w1i.dx(-1), w1i.dx(-2)),
                    //        dw2 = advection(c1, w2i.dx(-2), w2i.dx(-1), w2i, w2i.dx(+1), w2i.dx(+2)),
                    //        dw3 = advection(c2, w3i.dx(+2), w3i.dx(+1), w3i, w3i.dx(-1), w3i.dx(-2)),
                    //        dw4 = advection(c2, w4i.dx(-2), w4i.dx(-1), w4i, w4i.dx(+1), w4i.dx(+2)),

                    dw1 = advection(c1, w1i.dx(-2), w1i.dx(-1), w1i, w1i.dx(+1), w1i.dx(+2)),
                    dw2 = advection(c1, w2i.dx(+2), w2i.dx(+1), w2i, w2i.dx(-1), w2i.dx(-2)),
                    dw3 = advection(c2, w3i.dx(-2), w3i.dx(-1), w3i, w3i.dx(+1), w3i.dx(+2)),
                    dw4 = advection(c2, w4i.dx(+2), w4i.dx(+1), w4i, w4i.dx(-1), w4i.dx(-2)),

                    // from_omega_x
                    vxi += cpi * (dw1 - dw2),
                    vyi += csi * (dw3 - dw4),
                    sxxi += lmi * (dw1 + dw2),
                    syyi += lai * (dw1 + dw2),
                    sxyi += mui * (dw3 + dw4);
        }+
        [&]() { // to_omega_y
                    w1i = + vyi * icpi + syyi * ilmi,
                    w2i = - vyi * icpi + syyi * ilmi,
                    w3i = + vxi * icsi + sxyi * imui,
                    w4i = - vxi * icsi + sxyi * imui;
        }+
        [&]() { // omega_y_solve
                    c1 = cpi * dt / dy,
                    c2 = csi * dt / dy,

                    dw1 = advection(c1, w1i.dy(-2), w1i.dy(-1), w1i, w1i.dy(+1), w1i.dy(+2)),
                    dw2 = advection(c1, w2i.dy(+2), w2i.dy(+1), w2i, w2i.dy(-1), w2i.dy(-2)),
                    dw3 = advection(c2, w3i.dy(-2), w3i.dy(-1), w3i, w3i.dy(+1), w3i.dy(+2)),
                    dw4 = advection(c2, w4i.dy(+2), w4i.dy(+1), w4i, w4i.dy(-1), w4i.dy(-2)),

                    //        dw1 = advection(c1, w1i.dy(+2), w1i.dy(+1), w1i, w1i.dy(-1), w1i.dy(-2)),
                    //        dw2 = advection(c1, w2i.dy(-2), w2i.dy(-1), w2i, w2i.dy(+1), w2i.dy(+2)),
                    //        dw3 = advection(c2, w3i.dy(+2), w3i.dy(+1), w3i, w3i.dy(-1), w3i.dy(-2)),
                    //        dw4 = advection(c2, w4i.dy(-2), w4i.dy(-1), w4i, w4i.dy(+1), w4i.dy(+2)),

                    // from_omega_y
                    vxi += csi * (dw3 - dw4),
                    vyi += cpi * (dw1 - dw2),
                    sxxi += lai * (dw1 + dw2),
                    syyi += lmi * (dw1 + dw2),
                    sxyi += mui * (dw3 + dw4);
        };

        std::cout << "\rStep: " + std::to_string(t+1) << std::flush;

        if ((t + 1) % save_rate == 0) {
            vxa.save(save_path + "/vx_" + withLeadingZeros(t + 1, 6) + ".bin");
            vya.save(save_path + "/vy_" + withLeadingZeros(t + 1, 6) + ".bin");
            sxxa.save(save_path + "/sxx_" + withLeadingZeros(t + 1, 6) + ".bin");
            syya.save(save_path + "/syy_" + withLeadingZeros(t + 1, 6) + ".bin");
            sxya.save(save_path + "/sxy_" + withLeadingZeros(t + 1, 6) + ".bin");
        }
    }

    std::cout << std::endl;

    return 0;
}
