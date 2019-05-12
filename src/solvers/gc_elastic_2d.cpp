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
inline Expr le_max(Expr a, Expr b) { return a > b ? a : b; }
inline Expr le_max(Expr a, Expr b, Expr c) { return le_max(a, le_max(b, c)); }
inline bool isNull(Expr a) { return -RealEps < a && a < RealEps; }

Expr limiter(Expr r) {
    return le_max(0.0, le_min(1.0, 2.0 * r), le_min(2.0, r));
}

Expr advection(Expr c, Expr ppu, Expr pu, Expr u, Expr nu, Expr nnu) {
//    return - c * (u - pu);

    auto r1 = u - pu;
    auto r2 = nu - u;
    //if (isNull(r2)) r2 = RealEps;
    auto r = r1 / r2;

    auto pr1 = pu - ppu;
    auto pr2 = u - pu;
    //if (isNull(pr2)) pr2 = RealEps;
    auto pr = pr1 / pr2;

    auto nf12 = u + limiter(r) / 2.0 * (1.0 - c) * (nu - u);
    auto pf12 = pu + limiter(pr) / 2.0 * (1.0 - c) * (u - pu);

    return - c * (nf12 - pf12);
}

int main(int argc, char** argv) {
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
    Array cp_v(size), cs_v(size), rho_v(size), la_v(size), mu_v(size);
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

    {
        auto p = Array(nx, ny).load(conf["initial_pressure"].as<std::string>());
        sxx = Var(p);
        syy = Var(p);
    }

    {
        auto cpf = Array(nx, ny).load(conf["cp"].as<std::string>());
        cp = Var(cpf);
    }

    {
        auto csf = Array(nx, ny).load(conf["cs"].as<std::string>());
        cs = Var(csf);
    }

    {
        auto rhof = Array(nx, ny).load(conf["rho"].as<std::string>());
        rho = Var(rhof);
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

        mu = rho * cs * cs,
        la = rho * cp * cp - 2 * mu;
    }

    for (Int t = 0; t < time_steps; t++) {
        Temp c1, c2;
        Temp dw1, dw2, dw3, dw4;

        // to_omega_x
        w1 = + vx / (2 * cp) + sxx / (2 * (la + 2 * mu)),
        w2 = - vx / (2 * cp) + sxx / (2 * (la + 2 * mu)),
        w3 = + vy / (2 * cs) + sxy / (2 * mu),
        w4 = - vy / (2 * cs) + sxy / (2 * mu);

        // omega_x_solve
        c1 = cp * dt / dx,
        c2 = cs * dt / dx,
        dw1 = advection(c1, w1.dx(+2), w1.dx(+1), w1, w1.dx(-1), w1.dx(-2)),
        dw2 = advection(c1, w2.dx(-2), w2.dx(-1), w2, w2.dx(+1), w2.dx(+2)),
        dw3 = advection(c2, w3.dx(+2), w3.dx(+1), w3, w3.dx(-1), w3.dx(-2)),
        dw4 = advection(c2, w4.dx(-2), w4.dx(-1), w4, w4.dx(+1), w4.dx(+2)),

        // from_omega_x
        vx += cp * (dw1 - dw2),
        vy += cs * (dw3 - dw4),
        sxx += (la + 2 * mu) * (dw1 + dw2),
        syy += la * (dw1 + dw2),
        sxy += mu * (dw3 + dw4);

        // to_omega_y
        w1 = + vy / (2 * cp) + syy / (2 * (la + 2 * mu)),
        w2 = - vy / (2 * cp) + syy / (2 * (la + 2 * mu)),
        w3 = + vx / (2 * cs) + sxy / (2 * mu),
        w4 = - vx / (2 * cs) + sxy / (2 * mu);

        // omega_y_solve
        c1 = cp * dt / dy,
        c2 = cs * dt / dy,
        dw1 = advection(c1, w1.dy(+2), w1.dy(+1), w1, w1.dy(-1), w1.dy(-2)),
        dw2 = advection(c1, w2.dy(-2), w2.dy(-1), w2, w2.dy(+1), w2.dy(+2)),
        dw3 = advection(c2, w3.dy(+2), w3.dy(+1), w3, w3.dy(-1), w3.dy(-2)),
        dw4 = advection(c2, w4.dy(-2), w4.dy(-1), w4, w4.dy(+1), w4.dy(+2)),

        // from_omega_y
        vx += cs * (dw3 - dw4),
        vy += cp * (dw1 - dw2),
        sxx += la * (dw1 + dw2),
        syy += (la + 2 * mu) * (dw1 + dw2),
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
