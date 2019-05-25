#include "solverhelper.h"
#include "hrdata/hrdata.h"

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

        DefineLoop([time_steps, size](int Nx, int Ny) {
            for(Int j = 0; j < size.y; j++) {
                for(Int i = 0; i < size.x; i++) {
                    LoopSequence(0).init(i, j);
                    LoopSequence(1).init(i, j);
                    LoopSequence(2).init(i, j);
                    LoopSequence(3).init(i, j);
                }
            }
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
