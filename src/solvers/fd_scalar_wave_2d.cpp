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

    Array p_v(size), np_v(size), c_v(size);


    auto c = Var(c_v, inside);

    {
        auto pp = Array(nx, ny).load(conf["initial_pressure"].as<std::string>());
        Var(p_v, inside) = Var(pp);
    }

    {
        auto cc = Array(nx, ny).load(conf["cp"].as<std::string>());
        c = Var(cc);
    }


    // extend material values on the boundaries
    extend_ghost(Var(c_v), gh, gnx, gny);

    for (Int t = 0; t < time_steps; t++) {
        auto p = Var(t % 2 == 0 ? p_v : np_v, inside);
        auto np = Var(t % 2 == 0 ? np_v : p_v, inside);

        np = 2 * p - np + c * c * dt * dt * ((p.dx(-1) - 2 * p + p.dx(+1)) / (dx * dx) +
                                             (p.dy(-1) - 2 * p + p.dy(+1)) / (dy * dx));

        std::cout << "\rStep: " + std::to_string(t+1) << std::flush;

        if ((t + 1) % save_rate == 0) {
            p_v.save(save_path + "/sxx_" + withLeadingZeros(t + 1, 6) + ".bin");
        }
    }

    std::cout << std::endl;

    return 0;
}
