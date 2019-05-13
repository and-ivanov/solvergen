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

inline Expr le_min(Expr a, Expr b) { return a > b ? b : a; }
inline Expr le_max(Expr a, Expr b) { return a > b ? a : b;  }
inline Expr le_max(Expr a, Expr b, Expr c) { return le_max(a, le_max(b, c)); }
inline bool isNull(Expr a) { return -RealEps < a && a < RealEps; }

inline Expr limiter(Expr r) {
    return le_max(0.0, le_min(1.0, 2.0 * r), le_min(2.0, r));
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

    Array vx_v(size), vy_v(size), sxx_v(size), syy_v(size), sxy_v(size);


    Array cp_v(size), cs_v(size), rho_v(size), la_v(size), mu_v(size), lm_v(size);
    Array icp_v(size), ics_v(size), ila_v(size), imu_v(size), ilm_v(size);

    Array w1_v(size), w2_v(size), w3_v(size), w4_v(size);

    auto w1 = Var(w1_v, inside);
    auto w2 = Var(w2_v, inside);
    auto w3 = Var(w3_v, inside);
    auto w4 = Var(w4_v, inside);

    auto vx = Var(vx_v, inside);
    auto vy = Var(vy_v, inside);
    auto sxx = Var(sxx_v, inside);
    auto syy = Var(syy_v, inside);
    auto sxy = Var(sxy_v, inside);

    auto cp = Var(cp_v, inside);
    auto cs = Var(cs_v, inside);
    auto rho = Var(rho_v, inside);
    auto la = Var(la_v, inside);
    auto mu = Var(mu_v, inside);
    auto lm = Var(lm_v, inside);

    auto icp = Var(icp_v, inside);
    auto ics = Var(ics_v, inside);
    auto imu = Var(imu_v, inside);
    auto ilm = Var(ilm_v, inside);

    {
        auto p_v = Array(nx, ny).load(conf["initial_pressure"].as<std::string>());
        Var p(p_v);
        sxx = p;
        syy = p;
    }

    {
        auto cpf_v = Array(nx, ny).load(conf["cp"].as<std::string>());
        Var cpf(cpf_v);
        cp = cpf;
    }

    {
        auto csf_v = Array(nx, ny).load(conf["cs"].as<std::string>());
        Var csf(csf_v);
        cs = csf;
    }

    {
        auto rhof_v = Array(nx, ny).load(conf["rho"].as<std::string>());
        Var rhof(rhof_v);
        rho = rhof;
    }

    // extend material values on the boundaries
    extend_ghost(Var(cp_v), gh, gnx, gny);
    extend_ghost(Var(cs_v), gh, gnx, gny);
    extend_ghost(Var(rho_v), gh, gnx, gny);

    {
        auto cp = Var(cp_v);
        auto cs = Var(cs_v);
        auto rho = Var(rho_v);
        auto la = Var(la_v);
        auto mu = Var(mu_v);
        auto lm = Var(lm_v);

        auto icp = Var(icp_v);
        auto ics = Var(ics_v);
        auto imu = Var(imu_v);
        auto ilm = Var(ilm_v);

        mu = rho * cs * cs,
        la = rho * cp * cp - 2 * mu;

        lm = la + 2. * mu;

        icp = 1. / (2. * cp);
        ics = 1. / (2. * cs);
        imu = 1. / (2. * mu);
        ilm = 1. / (2. * lm);
    }

    Var w1nnx = w1.dx(+2);
    Var w1nx = w1.dx(+1);
    Var w1px = w1.dx(-1);
    Var w1ppx = w1.dx(-2);

    Var w1nny = w1.dy(+2);
    Var w1ny = w1.dy(+1);
    Var w1py = w1.dy(-1);
    Var w1ppy = w1.dy(-2);

    Var w2nnx = w2.dx(+2);
    Var w2nx = w2.dx(+1);
    Var w2px = w2.dx(-1);
    Var w2ppx = w2.dx(-2);

    Var w2nny = w2.dy(+2);
    Var w2ny = w2.dy(+1);
    Var w2py = w2.dy(-1);
    Var w2ppy = w2.dy(-2);

    Var w3nnx = w3.dx(+2);
    Var w3nx = w3.dx(+1);
    Var w3px = w3.dx(-1);
    Var w3ppx = w3.dx(-2);

    Var w3nny = w3.dy(+2);
    Var w3ny = w3.dy(+1);
    Var w3py = w3.dy(-1);
    Var w3ppy = w3.dy(-2);

    Var w4nnx = w4.dx(+2);
    Var w4nx = w4.dx(+1);
    Var w4px = w4.dx(-1);
    Var w4ppx = w4.dx(-2);

    Var w4nny = w4.dy(+2);
    Var w4ny = w4.dy(+1);
    Var w4py = w4.dy(-1);
    Var w4ppy = w4.dy(-2);

    for (Int t = 0; t < time_steps; t++) {
        Temp c1, c2;
        Temp dw1, dw2, dw3, dw4;

        // to_omega_x
        w1 = + vx * icp + sxx * ilm,
        w2 = - vx * icp + sxx * ilm,
        w3 = + vy * ics + sxy * imu,
        w4 = - vy * ics + sxy * imu;

        // omega_x_solve
        c1 = cp * dt / dx,
        c2 = cs * dt / dx,
        dw1 = advection(c1, w1nnx, w1nx, w1, w1px, w1ppx),
        dw2 = advection(c1, w2ppx, w2px, w2, w2nx, w2nnx),
        dw3 = advection(c2, w3nnx, w3nx, w3, w3px, w3ppx),
        dw4 = advection(c2, w4ppx, w4px, w4, w4nx, w4nnx),

        // from_omega_x
        vx += cp * (dw1 - dw2),
        vy += cs * (dw3 - dw4),
        sxx += lm * (dw1 + dw2),
        syy += la * (dw1 + dw2),
        sxy += mu * (dw3 + dw4);

        // to_omega_y
        w1 = + vy * icp + syy * ilm,
        w2 = - vy * icp + syy * ilm,
        w3 = + vx * ics + sxy * imu,
        w4 = - vx * ics + sxy * imu;

        // omega_y_solve
        c1 = cp * dt / dy,
        c2 = cs * dt / dy,
        dw1 = advection(c1, w1nny, w1ny, w1, w1py, w1ppy),
        dw2 = advection(c1, w2ppy, w2py, w2, w2ny, w2nny),
        dw3 = advection(c2, w3nny, w3ny, w3, w3py, w3ppy),
        dw4 = advection(c2, w4ppy, w4py, w4, w4ny, w4nny),

        // from_omega_y
        vx += cs * (dw3 - dw4),
        vy += cp * (dw1 - dw2),
        sxx += la * (dw1 + dw2),
        syy += lm * (dw1 + dw2),
        sxy += mu * (dw3 + dw4);

        std::cout << "\rStep: " + std::to_string(t+1) << std::flush;

        if ((t + 1) % save_rate == 0) {
            vx_v.save(save_path + "/vx_" + withLeadingZeros(t + 1, 6) + ".bin");
            vy_v.save(save_path + "/vy_" + withLeadingZeros(t + 1, 6) + ".bin");
            sxx_v.save(save_path + "/sxx_" + withLeadingZeros(t + 1, 6) + ".bin");
            syy_v.save(save_path + "/syy_" + withLeadingZeros(t + 1, 6) + ".bin");
            sxy_v.save(save_path + "/sxy_" + withLeadingZeros(t + 1, 6) + ".bin");
        }
    }

    std::cout << std::endl;

    return 0;
}
