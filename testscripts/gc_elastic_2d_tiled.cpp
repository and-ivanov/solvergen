#include "solverhelper.h"
#include "hrdata/hrdata.h"

#include <chrono>

void extend_ghost(Var a, Int gh, Int gnx, Int gny) {
    for (Int i = 0; i < gh; i++) {
        {
        Var t_a_x_gh_ghp1_ = a.x(gh,gh+1);
        Var t_a_x_i_ip1_ = a.x(i,i+1);
        Int _sizeX = t_a_x_gh_ghp1_.range().size().x;
        Int _sizeY = t_a_x_gh_ghp1_.range().size().y;
        assert(t_a_x_gh_ghp1_.range().size().x == _sizeX);
        assert(t_a_x_gh_ghp1_.range().size().y == _sizeY);
        assert(t_a_x_i_ip1_.range().size().x == _sizeX);
        assert(t_a_x_i_ip1_.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _a_x_gh_ghp1_;
        T _a_x_i_ip1_;
        // loads
        loadFromPtr(_a_x_gh_ghp1_, &t_a_x_gh_ghp1_.val(_i, _j));
        loadFromPtr(_a_x_i_ip1_, &t_a_x_i_ip1_.val(_i, _j));
        _a_x_i_ip1_ = + _a_x_gh_ghp1_;
        storeToPtr(&t_a_x_i_ip1_.val(_i, _j),_a_x_i_ip1_);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
        {
        Var t_a_x_gnxmghm1_gnxmgh_ = a.x(gnx-gh-1,gnx-gh);
        Var t_a_x_gnxmim1_gnxmi_ = a.x(gnx-i-1,gnx-i);
        Int _sizeX = t_a_x_gnxmghm1_gnxmgh_.range().size().x;
        Int _sizeY = t_a_x_gnxmghm1_gnxmgh_.range().size().y;
        assert(t_a_x_gnxmghm1_gnxmgh_.range().size().x == _sizeX);
        assert(t_a_x_gnxmghm1_gnxmgh_.range().size().y == _sizeY);
        assert(t_a_x_gnxmim1_gnxmi_.range().size().x == _sizeX);
        assert(t_a_x_gnxmim1_gnxmi_.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _a_x_gnxmghm1_gnxmgh_;
        T _a_x_gnxmim1_gnxmi_;
        // loads
        loadFromPtr(_a_x_gnxmghm1_gnxmgh_, &t_a_x_gnxmghm1_gnxmgh_.val(_i, _j));
        loadFromPtr(_a_x_gnxmim1_gnxmi_, &t_a_x_gnxmim1_gnxmi_.val(_i, _j));
        _a_x_gnxmim1_gnxmi_ = + _a_x_gnxmghm1_gnxmgh_;
        storeToPtr(&t_a_x_gnxmim1_gnxmi_.val(_i, _j),_a_x_gnxmim1_gnxmi_);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
    }
    for (Int j = 0; j < gh; j++) {
        {
        Var t_a_y_gh_ghp1_ = a.y(gh,gh+1);
        Var t_a_y_j_jp1_ = a.y(j,j+1);
        Int _sizeX = t_a_y_gh_ghp1_.range().size().x;
        Int _sizeY = t_a_y_gh_ghp1_.range().size().y;
        assert(t_a_y_gh_ghp1_.range().size().x == _sizeX);
        assert(t_a_y_gh_ghp1_.range().size().y == _sizeY);
        assert(t_a_y_j_jp1_.range().size().x == _sizeX);
        assert(t_a_y_j_jp1_.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _a_y_gh_ghp1_;
        T _a_y_j_jp1_;
        // loads
        loadFromPtr(_a_y_gh_ghp1_, &t_a_y_gh_ghp1_.val(_i, _j));
        loadFromPtr(_a_y_j_jp1_, &t_a_y_j_jp1_.val(_i, _j));
        _a_y_j_jp1_ = + _a_y_gh_ghp1_;
        storeToPtr(&t_a_y_j_jp1_.val(_i, _j),_a_y_j_jp1_);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
        {
        Var t_a_y_gnymghm1_gnymgh_ = a.y(gny-gh-1,gny-gh);
        Var t_a_y_gnymjm1_gnymj_ = a.y(gny-j-1,gny-j);
        Int _sizeX = t_a_y_gnymghm1_gnymgh_.range().size().x;
        Int _sizeY = t_a_y_gnymghm1_gnymgh_.range().size().y;
        assert(t_a_y_gnymghm1_gnymgh_.range().size().x == _sizeX);
        assert(t_a_y_gnymghm1_gnymgh_.range().size().y == _sizeY);
        assert(t_a_y_gnymjm1_gnymj_.range().size().x == _sizeX);
        assert(t_a_y_gnymjm1_gnymj_.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _a_y_gnymghm1_gnymgh_;
        T _a_y_gnymjm1_gnymj_;
        // loads
        loadFromPtr(_a_y_gnymghm1_gnymgh_, &t_a_y_gnymghm1_gnymgh_.val(_i, _j));
        loadFromPtr(_a_y_gnymjm1_gnymj_, &t_a_y_gnymjm1_gnymj_.val(_i, _j));
        _a_y_gnymjm1_gnymj_ = + _a_y_gnymghm1_gnymgh_;
        storeToPtr(&t_a_y_gnymjm1_gnymj_.val(_i, _j),_a_y_gnymjm1_gnymj_);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
    }
}

template <typename E0>
__attribute__((always_inline,flatten))
inline bool isNull(E0 a) { return -RealEps < a && a < RealEps; }

template <typename E0>
__attribute__((always_inline,flatten))
inline auto limiter(E0 r) {
    return max(0.0, max(min(1.0, 2.0 * r), min(2.0, r)));
}

template <typename E0, typename E1, typename E2, typename E3, typename E4, typename E5>
__attribute__((always_inline,flatten))
inline auto advection(E0 c, E1 ppu, E2 pu, E3 u, E4 nu, E5 nnu)
{
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

    Int _sizeX = gnx;
Int _sizeY = gny;


    // last letter a denotes entire Array
    // last letter i of Var denotes inner area of Array without ghost

    Array vxa(_sizeX, _sizeY), vya(_sizeX, _sizeY), sxxa(_sizeX, _sizeY), syya(_sizeX, _sizeY), sxya(_sizeX, _sizeY);


    Array cpa(_sizeX, _sizeY), csa(_sizeX, _sizeY), rhoa(_sizeX, _sizeY), laa(_sizeX, _sizeY), mua(_sizeX, _sizeY), lma(_sizeX, _sizeY);
    Array icpa(_sizeX, _sizeY), icsa(_sizeX, _sizeY), ilaa(_sizeX, _sizeY), imua(_sizeX, _sizeY), ilma(_sizeX, _sizeY);

    Array w1a(_sizeX, _sizeY), w2a(_sizeX, _sizeY), w3a(_sizeX, _sizeY), w4a(_sizeX, _sizeY);

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
        {
        Var t_pi = pi;
        Var t_sxxi = sxxi;
        Int _sizeX = t_pi.range().size().x;
        Int _sizeY = t_pi.range().size().y;
        assert(t_pi.range().size().x == _sizeX);
        assert(t_pi.range().size().y == _sizeY);
        assert(t_sxxi.range().size().x == _sizeX);
        assert(t_sxxi.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _pi;
        T _sxxi;
        // loads
        loadFromPtr(_pi, &t_pi.val(_i, _j));
        loadFromPtr(_sxxi, &t_sxxi.val(_i, _j));
        _sxxi = _pi;
        storeToPtr(&t_sxxi.val(_i, _j),_sxxi);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
        {
        Var t_pi = pi;
        Var t_syyi = syyi;
        Int _sizeX = t_pi.range().size().x;
        Int _sizeY = t_pi.range().size().y;
        assert(t_pi.range().size().x == _sizeX);
        assert(t_pi.range().size().y == _sizeY);
        assert(t_syyi.range().size().x == _sizeX);
        assert(t_syyi.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _pi;
        T _syyi;
        // loads
        loadFromPtr(_pi, &t_pi.val(_i, _j));
        loadFromPtr(_syyi, &t_syyi.val(_i, _j));
        _syyi = _pi;
        storeToPtr(&t_syyi.val(_i, _j),_syyi);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
    }

    {
        auto cpf = Array(nx, ny).load(conf["cp"].as<std::string>());
        {
        Var t_cpf = cpf;
        Var t_cpi = cpi;
        Int _sizeX = t_cpf.range().size().x;
        Int _sizeY = t_cpf.range().size().y;
        assert(t_cpf.range().size().x == _sizeX);
        assert(t_cpf.range().size().y == _sizeY);
        assert(t_cpi.range().size().x == _sizeX);
        assert(t_cpi.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _cpf;
        T _cpi;
        // loads
        loadFromPtr(_cpf, &t_cpf.val(_i, _j));
        loadFromPtr(_cpi, &t_cpi.val(_i, _j));
        _cpi = _cpf;
        storeToPtr(&t_cpi.val(_i, _j),_cpi);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
    }

    {
        auto csf = Array(nx, ny).load(conf["cs"].as<std::string>());
        {
        Var t_csf = csf;
        Var t_csi = csi;
        Int _sizeX = t_csf.range().size().x;
        Int _sizeY = t_csf.range().size().y;
        assert(t_csf.range().size().x == _sizeX);
        assert(t_csf.range().size().y == _sizeY);
        assert(t_csi.range().size().x == _sizeX);
        assert(t_csi.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _csf;
        T _csi;
        // loads
        loadFromPtr(_csf, &t_csf.val(_i, _j));
        loadFromPtr(_csi, &t_csi.val(_i, _j));
        _csi = _csf;
        storeToPtr(&t_csi.val(_i, _j),_csi);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
    }

    {
        auto rhof = Array(nx, ny).load(conf["rho"].as<std::string>());
        {
        Var t_rhof = rhof;
        Var t_rhoi = rhoi;
        Int _sizeX = t_rhof.range().size().x;
        Int _sizeY = t_rhof.range().size().y;
        assert(t_rhof.range().size().x == _sizeX);
        assert(t_rhof.range().size().y == _sizeY);
        assert(t_rhoi.range().size().x == _sizeX);
        assert(t_rhoi.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _rhof;
        T _rhoi;
        // loads
        loadFromPtr(_rhof, &t_rhof.val(_i, _j));
        loadFromPtr(_rhoi, &t_rhoi.val(_i, _j));
        _rhoi = _rhof;
        storeToPtr(&t_rhoi.val(_i, _j),_rhoi);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
    }

    // extend material values on the boundaries
    extend_ghost(cpa, gh, gnx, gny);
    extend_ghost(csa, gh, gnx, gny);
    extend_ghost(rhoa, gh, gnx, gny);

    {
        {
        Var t_cpa = cpa;
        Var t_csa = csa;
        Var t_laa = laa;
        Var t_mua = mua;
        Var t_rhoa = rhoa;
        Int _sizeX = t_cpa.range().size().x;
        Int _sizeY = t_cpa.range().size().y;
        assert(t_cpa.range().size().x == _sizeX);
        assert(t_cpa.range().size().y == _sizeY);
        assert(t_csa.range().size().x == _sizeX);
        assert(t_csa.range().size().y == _sizeY);
        assert(t_laa.range().size().x == _sizeX);
        assert(t_laa.range().size().y == _sizeY);
        assert(t_mua.range().size().x == _sizeX);
        assert(t_mua.range().size().y == _sizeY);
        assert(t_rhoa.range().size().x == _sizeX);
        assert(t_rhoa.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _cpa;
        T _csa;
        T _laa;
        T _mua;
        T _rhoa;
        // loads
        loadFromPtr(_cpa, &t_cpa.val(_i, _j));
        loadFromPtr(_csa, &t_csa.val(_i, _j));
        loadFromPtr(_laa, &t_laa.val(_i, _j));
        loadFromPtr(_mua, &t_mua.val(_i, _j));
        loadFromPtr(_rhoa, &t_rhoa.val(_i, _j));
        _mua = _rhoa * _csa * _csa;
        _laa = _rhoa * _cpa * _cpa - 2 * _mua;
        storeToPtr(&t_laa.val(_i, _j),_laa);
        storeToPtr(&t_mua.val(_i, _j),_mua);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        

        {
        Var t_laa = laa;
        Var t_lma = lma;
        Var t_mua = mua;
        Int _sizeX = t_laa.range().size().x;
        Int _sizeY = t_laa.range().size().y;
        assert(t_laa.range().size().x == _sizeX);
        assert(t_laa.range().size().y == _sizeY);
        assert(t_lma.range().size().x == _sizeX);
        assert(t_lma.range().size().y == _sizeY);
        assert(t_mua.range().size().x == _sizeX);
        assert(t_mua.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _laa;
        T _lma;
        T _mua;
        // loads
        loadFromPtr(_laa, &t_laa.val(_i, _j));
        loadFromPtr(_lma, &t_lma.val(_i, _j));
        loadFromPtr(_mua, &t_mua.val(_i, _j));
        _lma = _laa + 2. * _mua;
        storeToPtr(&t_lma.val(_i, _j),_lma);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        

        {
        Var t_cpa = cpa;
        Var t_icpa = icpa;
        Int _sizeX = t_cpa.range().size().x;
        Int _sizeY = t_cpa.range().size().y;
        assert(t_cpa.range().size().x == _sizeX);
        assert(t_cpa.range().size().y == _sizeY);
        assert(t_icpa.range().size().x == _sizeX);
        assert(t_icpa.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _cpa;
        T _icpa;
        // loads
        loadFromPtr(_cpa, &t_cpa.val(_i, _j));
        loadFromPtr(_icpa, &t_icpa.val(_i, _j));
        _icpa = 1. / (2. * _cpa);
        storeToPtr(&t_icpa.val(_i, _j),_icpa);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
        {
        Var t_csa = csa;
        Var t_icsa = icsa;
        Int _sizeX = t_csa.range().size().x;
        Int _sizeY = t_csa.range().size().y;
        assert(t_csa.range().size().x == _sizeX);
        assert(t_csa.range().size().y == _sizeY);
        assert(t_icsa.range().size().x == _sizeX);
        assert(t_icsa.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _csa;
        T _icsa;
        // loads
        loadFromPtr(_csa, &t_csa.val(_i, _j));
        loadFromPtr(_icsa, &t_icsa.val(_i, _j));
        _icsa = 1. / (2. * _csa);
        storeToPtr(&t_icsa.val(_i, _j),_icsa);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
        {
        Var t_imua = imua;
        Var t_mua = mua;
        Int _sizeX = t_imua.range().size().x;
        Int _sizeY = t_imua.range().size().y;
        assert(t_imua.range().size().x == _sizeX);
        assert(t_imua.range().size().y == _sizeY);
        assert(t_mua.range().size().x == _sizeX);
        assert(t_mua.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _imua;
        T _mua;
        // loads
        loadFromPtr(_imua, &t_imua.val(_i, _j));
        loadFromPtr(_mua, &t_mua.val(_i, _j));
        _imua = 1. / (2. * _mua);
        storeToPtr(&t_imua.val(_i, _j),_imua);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
        {
        Var t_ilma = ilma;
        Var t_lma = lma;
        Int _sizeX = t_ilma.range().size().x;
        Int _sizeY = t_ilma.range().size().y;
        assert(t_ilma.range().size().x == _sizeX);
        assert(t_ilma.range().size().y == _sizeY);
        assert(t_lma.range().size().x == _sizeX);
        assert(t_lma.range().size().y == _sizeY);
        auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {
        using T = decltype(type);
        // temporaries
        // variables
        T _ilma;
        T _lma;
        // loads
        loadFromPtr(_ilma, &t_ilma.val(_i, _j));
        loadFromPtr(_lma, &t_lma.val(_i, _j));
        _ilma = 1. / (2. * _lma);
        storeToPtr(&t_ilma.val(_i, _j),_ilma);
        };
        for (Int _j = 0; _j < _sizeY; _j++) {
        Int _i = 0;
        for (; _i < _sizeX / 4 * 4; _i+=4) {
        _evaluate(_i, _j, Real4());
        }
        for (; _i < _sizeX; _i++) {
        _evaluate(_i, _j, Real());
        }
        }
        }
        
    }

    //for (Int t = 0; t < time_steps; t++) {
    {
        // to_omega_x
        Var t_icpi = icpi;
        Var t_icsi = icsi;
        Var t_ilmi = ilmi;
        Var t_imui = imui;
        Var t_sxxi = sxxi;
        Var t_sxyi = sxyi;
        Var t_vxi = vxi;
        Var t_vyi = vyi;
        Var t_w1i = w1i;
        Var t_w2i = w2i;
        Var t_w3i = w3i;
        Var t_w4i = w4i;
        Int _sizeX = t_icpi.range().size().x;
        Int _sizeY = t_icpi.range().size().y;
        assert(t_icpi.range().size().x == _sizeX);
        assert(t_icpi.range().size().y == _sizeY);
        assert(t_icsi.range().size().x == _sizeX);
        assert(t_icsi.range().size().y == _sizeY);
        assert(t_ilmi.range().size().x == _sizeX);
        assert(t_ilmi.range().size().y == _sizeY);
        assert(t_imui.range().size().x == _sizeX);
        assert(t_imui.range().size().y == _sizeY);
        assert(t_sxxi.range().size().x == _sizeX);
        assert(t_sxxi.range().size().y == _sizeY);
        assert(t_sxyi.range().size().x == _sizeX);
        assert(t_sxyi.range().size().y == _sizeY);
        assert(t_vxi.range().size().x == _sizeX);
        assert(t_vxi.range().size().y == _sizeY);
        assert(t_vyi.range().size().x == _sizeX);
        assert(t_vyi.range().size().y == _sizeY);
        assert(t_w1i.range().size().x == _sizeX);
        assert(t_w1i.range().size().y == _sizeY);
        assert(t_w2i.range().size().x == _sizeX);
        assert(t_w2i.range().size().y == _sizeY);
        assert(t_w3i.range().size().x == _sizeX);
        assert(t_w3i.range().size().y == _sizeY);
        assert(t_w4i.range().size().x == _sizeX);
        assert(t_w4i.range().size().y == _sizeY);
        auto _evaluateX1 =
                [&]
                (Int _i, Int _j, auto type)
                __attribute__((always_inline,flatten))
        {
        using T = decltype(type);
        // temporaries
        // variables
        T _icpi;
        T _icsi;
        T _ilmi;
        T _imui;
        T _sxxi;
        T _sxyi;
        T _vxi;
        T _vyi;
        T _w1i;
        T _w2i;
        T _w3i;
        T _w4i;
        // loads
        loadFromPtr(_icpi, &t_icpi.val(_i, _j));
        loadFromPtr(_icsi, &t_icsi.val(_i, _j));
        loadFromPtr(_ilmi, &t_ilmi.val(_i, _j));
        loadFromPtr(_imui, &t_imui.val(_i, _j));
        loadFromPtr(_sxxi, &t_sxxi.val(_i, _j));
        loadFromPtr(_sxyi, &t_sxyi.val(_i, _j));
        loadFromPtr(_vxi, &t_vxi.val(_i, _j));
        loadFromPtr(_vyi, &t_vyi.val(_i, _j));
        loadFromPtr(_w1i, &t_w1i.val(_i, _j));
        loadFromPtr(_w2i, &t_w2i.val(_i, _j));
        loadFromPtr(_w3i, &t_w3i.val(_i, _j));
        loadFromPtr(_w4i, &t_w4i.val(_i, _j));
        _w1i = + _vxi * _icpi + _sxxi * _ilmi;
        _w2i = - _vxi * _icpi + _sxxi * _ilmi;
        _w3i = + _vyi * _icsi + _sxyi * _imui;
        _w4i = - _vyi * _icsi + _sxyi * _imui;
        storeToPtr(&t_w1i.val(_i, _j),_w1i);
        storeToPtr(&t_w2i.val(_i, _j),_w2i);
        storeToPtr(&t_w3i.val(_i, _j),_w3i);
        storeToPtr(&t_w4i.val(_i, _j),_w4i);
        };
//        for (Int _j = 0; _j < _sizeY; _j++) {
//        Int _i = 0;
//        for (; _i < _sizeX / 4 * 4; _i+=4) {
//        _evaluateX1(_i, _j, Real4());
//        }
//        for (; _i < _sizeX; _i++) {
//        _evaluateX1(_i, _j, Real());
//        }
//        }
        

        // omega_x_solve
        Var t_cpi = cpi;
        Var t_csi = csi;
        Var t_lai = lai;
        Var t_lmi = lmi;
        Var t_mui = mui;
//        Var t_sxxi = sxxi;
//        Var t_sxyi = sxyi;
        Var t_syyi = syyi;
//        Var t_vxi = vxi;
//        Var t_vyi = vyi;
//        Var t_w1i = w1i;
        Var t_w1i_dx_p1_ = w1i.dx(+1);
        Var t_w1i_dx_p2_ = w1i.dx(+2);
        Var t_w1i_dx_m1_ = w1i.dx(-1);
        Var t_w1i_dx_m2_ = w1i.dx(-2);
//        Var t_w2i = w2i;
        Var t_w2i_dx_p1_ = w2i.dx(+1);
        Var t_w2i_dx_p2_ = w2i.dx(+2);
        Var t_w2i_dx_m1_ = w2i.dx(-1);
        Var t_w2i_dx_m2_ = w2i.dx(-2);
//        Var t_w3i = w3i;
        Var t_w3i_dx_p1_ = w3i.dx(+1);
        Var t_w3i_dx_p2_ = w3i.dx(+2);
        Var t_w3i_dx_m1_ = w3i.dx(-1);
        Var t_w3i_dx_m2_ = w3i.dx(-2);
//        Var t_w4i = w4i;
        Var t_w4i_dx_p1_ = w4i.dx(+1);
        Var t_w4i_dx_p2_ = w4i.dx(+2);
        Var t_w4i_dx_m1_ = w4i.dx(-1);
        Var t_w4i_dx_m2_ = w4i.dx(-2);
//        Int _sizeX = t_cpi.range().size().x;
//        Int _sizeY = t_cpi.range().size().y;
        assert(t_cpi.range().size().x == _sizeX);
        assert(t_cpi.range().size().y == _sizeY);
        assert(t_csi.range().size().x == _sizeX);
        assert(t_csi.range().size().y == _sizeY);
        assert(t_lai.range().size().x == _sizeX);
        assert(t_lai.range().size().y == _sizeY);
        assert(t_lmi.range().size().x == _sizeX);
        assert(t_lmi.range().size().y == _sizeY);
        assert(t_mui.range().size().x == _sizeX);
        assert(t_mui.range().size().y == _sizeY);
        assert(t_sxxi.range().size().x == _sizeX);
        assert(t_sxxi.range().size().y == _sizeY);
        assert(t_sxyi.range().size().x == _sizeX);
        assert(t_sxyi.range().size().y == _sizeY);
        assert(t_syyi.range().size().x == _sizeX);
        assert(t_syyi.range().size().y == _sizeY);
        assert(t_vxi.range().size().x == _sizeX);
        assert(t_vxi.range().size().y == _sizeY);
        assert(t_vyi.range().size().x == _sizeX);
        assert(t_vyi.range().size().y == _sizeY);
        assert(t_w1i.range().size().x == _sizeX);
        assert(t_w1i.range().size().y == _sizeY);
        assert(t_w1i_dx_p1_.range().size().x == _sizeX);
        assert(t_w1i_dx_p1_.range().size().y == _sizeY);
        assert(t_w1i_dx_p2_.range().size().x == _sizeX);
        assert(t_w1i_dx_p2_.range().size().y == _sizeY);
        assert(t_w1i_dx_m1_.range().size().x == _sizeX);
        assert(t_w1i_dx_m1_.range().size().y == _sizeY);
        assert(t_w1i_dx_m2_.range().size().x == _sizeX);
        assert(t_w1i_dx_m2_.range().size().y == _sizeY);
        assert(t_w2i.range().size().x == _sizeX);
        assert(t_w2i.range().size().y == _sizeY);
        assert(t_w2i_dx_p1_.range().size().x == _sizeX);
        assert(t_w2i_dx_p1_.range().size().y == _sizeY);
        assert(t_w2i_dx_p2_.range().size().x == _sizeX);
        assert(t_w2i_dx_p2_.range().size().y == _sizeY);
        assert(t_w2i_dx_m1_.range().size().x == _sizeX);
        assert(t_w2i_dx_m1_.range().size().y == _sizeY);
        assert(t_w2i_dx_m2_.range().size().x == _sizeX);
        assert(t_w2i_dx_m2_.range().size().y == _sizeY);
        assert(t_w3i.range().size().x == _sizeX);
        assert(t_w3i.range().size().y == _sizeY);
        assert(t_w3i_dx_p1_.range().size().x == _sizeX);
        assert(t_w3i_dx_p1_.range().size().y == _sizeY);
        assert(t_w3i_dx_p2_.range().size().x == _sizeX);
        assert(t_w3i_dx_p2_.range().size().y == _sizeY);
        assert(t_w3i_dx_m1_.range().size().x == _sizeX);
        assert(t_w3i_dx_m1_.range().size().y == _sizeY);
        assert(t_w3i_dx_m2_.range().size().x == _sizeX);
        assert(t_w3i_dx_m2_.range().size().y == _sizeY);
        assert(t_w4i.range().size().x == _sizeX);
        assert(t_w4i.range().size().y == _sizeY);
        assert(t_w4i_dx_p1_.range().size().x == _sizeX);
        assert(t_w4i_dx_p1_.range().size().y == _sizeY);
        assert(t_w4i_dx_p2_.range().size().x == _sizeX);
        assert(t_w4i_dx_p2_.range().size().y == _sizeY);
        assert(t_w4i_dx_m1_.range().size().x == _sizeX);
        assert(t_w4i_dx_m1_.range().size().y == _sizeY);
        assert(t_w4i_dx_m2_.range().size().x == _sizeX);
        assert(t_w4i_dx_m2_.range().size().y == _sizeY);
        auto _evaluateX2 =
                [&] (Int _i, Int _j, auto type)
                __attribute__((always_inline,flatten))
        {
        using T = decltype(type);
        // temporaries
        T _c1;
        T _c2;
        T _dw1;
        T _dw2;
        T _dw3;
        T _dw4;
        // variables
        T _cpi;
        T _csi;
        T _lai;
        T _lmi;
        T _mui;
        T _sxxi;
        T _sxyi;
        T _syyi;
        T _vxi;
        T _vyi;
        T _w1i;
        T _w1i_dx_p1_;
        T _w1i_dx_p2_;
        T _w1i_dx_m1_;
        T _w1i_dx_m2_;
        T _w2i;
        T _w2i_dx_p1_;
        T _w2i_dx_p2_;
        T _w2i_dx_m1_;
        T _w2i_dx_m2_;
        T _w3i;
        T _w3i_dx_p1_;
        T _w3i_dx_p2_;
        T _w3i_dx_m1_;
        T _w3i_dx_m2_;
        T _w4i;
        T _w4i_dx_p1_;
        T _w4i_dx_p2_;
        T _w4i_dx_m1_;
        T _w4i_dx_m2_;
        // loads
        loadFromPtr(_cpi, &t_cpi.val(_i, _j));
        loadFromPtr(_csi, &t_csi.val(_i, _j));
        loadFromPtr(_lai, &t_lai.val(_i, _j));
        loadFromPtr(_lmi, &t_lmi.val(_i, _j));
        loadFromPtr(_mui, &t_mui.val(_i, _j));
        loadFromPtr(_sxxi, &t_sxxi.val(_i, _j));
        loadFromPtr(_sxyi, &t_sxyi.val(_i, _j));
        loadFromPtr(_syyi, &t_syyi.val(_i, _j));
        loadFromPtr(_vxi, &t_vxi.val(_i, _j));
        loadFromPtr(_vyi, &t_vyi.val(_i, _j));
        loadFromPtr(_w1i, &t_w1i.val(_i, _j));
        loadFromPtr(_w1i_dx_p1_, &t_w1i_dx_p1_.val(_i, _j));
        loadFromPtr(_w1i_dx_p2_, &t_w1i_dx_p2_.val(_i, _j));
        loadFromPtr(_w1i_dx_m1_, &t_w1i_dx_m1_.val(_i, _j));
        loadFromPtr(_w1i_dx_m2_, &t_w1i_dx_m2_.val(_i, _j));
        loadFromPtr(_w2i, &t_w2i.val(_i, _j));
        loadFromPtr(_w2i_dx_p1_, &t_w2i_dx_p1_.val(_i, _j));
        loadFromPtr(_w2i_dx_p2_, &t_w2i_dx_p2_.val(_i, _j));
        loadFromPtr(_w2i_dx_m1_, &t_w2i_dx_m1_.val(_i, _j));
        loadFromPtr(_w2i_dx_m2_, &t_w2i_dx_m2_.val(_i, _j));
        loadFromPtr(_w3i, &t_w3i.val(_i, _j));
        loadFromPtr(_w3i_dx_p1_, &t_w3i_dx_p1_.val(_i, _j));
        loadFromPtr(_w3i_dx_p2_, &t_w3i_dx_p2_.val(_i, _j));
        loadFromPtr(_w3i_dx_m1_, &t_w3i_dx_m1_.val(_i, _j));
        loadFromPtr(_w3i_dx_m2_, &t_w3i_dx_m2_.val(_i, _j));
        loadFromPtr(_w4i, &t_w4i.val(_i, _j));
        loadFromPtr(_w4i_dx_p1_, &t_w4i_dx_p1_.val(_i, _j));
        loadFromPtr(_w4i_dx_p2_, &t_w4i_dx_p2_.val(_i, _j));
        loadFromPtr(_w4i_dx_m1_, &t_w4i_dx_m1_.val(_i, _j));
        loadFromPtr(_w4i_dx_m2_, &t_w4i_dx_m2_.val(_i, _j));
        _c1 = _cpi * dt / dx;
        _c2 = _csi * dt / dx;

        _dw1 = advection(_c1, _w1i_dx_m2_, _w1i_dx_m1_, _w1i, _w1i_dx_p1_, _w1i_dx_p2_);
        _dw2 = advection(_c1, _w2i_dx_p2_, _w2i_dx_p1_, _w2i, _w2i_dx_m1_, _w2i_dx_m2_);
        _dw3 = advection(_c2, _w3i_dx_m2_, _w3i_dx_m1_, _w3i, _w3i_dx_p1_, _w3i_dx_p2_);
        _dw4 = advection(_c2, _w4i_dx_p2_, _w4i_dx_p1_, _w4i, _w4i_dx_m1_, _w4i_dx_m2_);

        // from_omega_x
        _vxi += _cpi * (_dw1 - _dw2);
        _vyi += _csi * (_dw3 - _dw4);
        _sxxi += _lmi * (_dw1 + _dw2);
        _syyi += _lai * (_dw1 + _dw2);
        _sxyi += _mui * (_dw3 + _dw4);
        storeToPtr(&t_sxxi.val(_i, _j),_sxxi);
        storeToPtr(&t_sxyi.val(_i, _j),_sxyi);
        storeToPtr(&t_syyi.val(_i, _j),_syyi);
        storeToPtr(&t_vxi.val(_i, _j),_vxi);
        storeToPtr(&t_vyi.val(_i, _j),_vyi);
        };
//        for (Int _j = 0; _j < _sizeY; _j++) {
//        Int _i = 0;
//        for (; _i < _sizeX / 4 * 4; _i+=4) {
//        _evaluateX2(_i, _j, Real4());
//        }
//        for (; _i < _sizeX; _i++) {
//        _evaluateX2(_i, _j, Real());
//        }
//        }
        

        // to_omega_y
//        Var t_icpi = icpi;
//        Var t_icsi = icsi;
//        Var t_ilmi = ilmi;
//        Var t_imui = imui;
//        Var t_sxyi = sxyi;
//        Var t_syyi = syyi;
//        Var t_vxi = vxi;
//        Var t_vyi = vyi;
//        Var t_w1i = w1i;
//        Var t_w2i = w2i;
//        Var t_w3i = w3i;
//        Var t_w4i = w4i;
//        Int _sizeX = t_icpi.range().size().x;
//        Int _sizeY = t_icpi.range().size().y;
        assert(t_icpi.range().size().x == _sizeX);
        assert(t_icpi.range().size().y == _sizeY);
        assert(t_icsi.range().size().x == _sizeX);
        assert(t_icsi.range().size().y == _sizeY);
        assert(t_ilmi.range().size().x == _sizeX);
        assert(t_ilmi.range().size().y == _sizeY);
        assert(t_imui.range().size().x == _sizeX);
        assert(t_imui.range().size().y == _sizeY);
        assert(t_sxyi.range().size().x == _sizeX);
        assert(t_sxyi.range().size().y == _sizeY);
        assert(t_syyi.range().size().x == _sizeX);
        assert(t_syyi.range().size().y == _sizeY);
        assert(t_vxi.range().size().x == _sizeX);
        assert(t_vxi.range().size().y == _sizeY);
        assert(t_vyi.range().size().x == _sizeX);
        assert(t_vyi.range().size().y == _sizeY);
        assert(t_w1i.range().size().x == _sizeX);
        assert(t_w1i.range().size().y == _sizeY);
        assert(t_w2i.range().size().x == _sizeX);
        assert(t_w2i.range().size().y == _sizeY);
        assert(t_w3i.range().size().x == _sizeX);
        assert(t_w3i.range().size().y == _sizeY);
        assert(t_w4i.range().size().x == _sizeX);
        assert(t_w4i.range().size().y == _sizeY);
        auto _evaluateY1 = [&]
                (Int _i, Int _j, auto type)
                __attribute__((always_inline,flatten))
        {
        using T = decltype(type);
        // temporaries
        // variables
        T _icpi;
        T _icsi;
        T _ilmi;
        T _imui;
        T _sxyi;
        T _syyi;
        T _vxi;
        T _vyi;
        T _w1i;
        T _w2i;
        T _w3i;
        T _w4i;
        // loads
        loadFromPtr(_icpi, &t_icpi.val(_i, _j));
        loadFromPtr(_icsi, &t_icsi.val(_i, _j));
        loadFromPtr(_ilmi, &t_ilmi.val(_i, _j));
        loadFromPtr(_imui, &t_imui.val(_i, _j));
        loadFromPtr(_sxyi, &t_sxyi.val(_i, _j));
        loadFromPtr(_syyi, &t_syyi.val(_i, _j));
        loadFromPtr(_vxi, &t_vxi.val(_i, _j));
        loadFromPtr(_vyi, &t_vyi.val(_i, _j));
        loadFromPtr(_w1i, &t_w1i.val(_i, _j));
        loadFromPtr(_w2i, &t_w2i.val(_i, _j));
        loadFromPtr(_w3i, &t_w3i.val(_i, _j));
        loadFromPtr(_w4i, &t_w4i.val(_i, _j));
        _w1i = + _vyi * _icpi + _syyi * _ilmi;
        _w2i = - _vyi * _icpi + _syyi * _ilmi;
        _w3i = + _vxi * _icsi + _sxyi * _imui;
        _w4i = - _vxi * _icsi + _sxyi * _imui;
        storeToPtr(&t_w1i.val(_i, _j),_w1i);
        storeToPtr(&t_w2i.val(_i, _j),_w2i);
        storeToPtr(&t_w3i.val(_i, _j),_w3i);
        storeToPtr(&t_w4i.val(_i, _j),_w4i);
        };
//        for (Int _j = 0; _j < _sizeY; _j++) {
//        Int _i = 0;
//        for (; _i < _sizeX / 4 * 4; _i+=4) {
//        _evaluateY1(_i, _j, Real4());
//        }
//        for (; _i < _sizeX; _i++) {
//        _evaluateY1(_i, _j, Real());
//        }
//        }
        

        // omega_y_solve
//        Var t_cpi = cpi;
//        Var t_csi = csi;
//        Var t_lai = lai;
//        Var t_lmi = lmi;
//        Var t_mui = mui;
//        Var t_sxxi = sxxi;
//        Var t_sxyi = sxyi;
//        Var t_syyi = syyi;
//        Var t_vxi = vxi;
//        Var t_vyi = vyi;
//        Var t_w1i = w1i;
        Var t_w1i_dy_p1_ = w1i.dy(+1);
        Var t_w1i_dy_p2_ = w1i.dy(+2);
        Var t_w1i_dy_m1_ = w1i.dy(-1);
        Var t_w1i_dy_m2_ = w1i.dy(-2);
//        Var t_w2i = w2i;
        Var t_w2i_dy_p1_ = w2i.dy(+1);
        Var t_w2i_dy_p2_ = w2i.dy(+2);
        Var t_w2i_dy_m1_ = w2i.dy(-1);
        Var t_w2i_dy_m2_ = w2i.dy(-2);
//        Var t_w3i = w3i;
        Var t_w3i_dy_p1_ = w3i.dy(+1);
        Var t_w3i_dy_p2_ = w3i.dy(+2);
        Var t_w3i_dy_m1_ = w3i.dy(-1);
        Var t_w3i_dy_m2_ = w3i.dy(-2);
//        Var t_w4i = w4i;
        Var t_w4i_dy_p1_ = w4i.dy(+1);
        Var t_w4i_dy_p2_ = w4i.dy(+2);
        Var t_w4i_dy_m1_ = w4i.dy(-1);
        Var t_w4i_dy_m2_ = w4i.dy(-2);
//        Int _sizeX = t_cpi.range().size().x;
//        Int _sizeY = t_cpi.range().size().y;
        assert(t_cpi.range().size().x == _sizeX);
        assert(t_cpi.range().size().y == _sizeY);
        assert(t_csi.range().size().x == _sizeX);
        assert(t_csi.range().size().y == _sizeY);
        assert(t_lai.range().size().x == _sizeX);
        assert(t_lai.range().size().y == _sizeY);
        assert(t_lmi.range().size().x == _sizeX);
        assert(t_lmi.range().size().y == _sizeY);
        assert(t_mui.range().size().x == _sizeX);
        assert(t_mui.range().size().y == _sizeY);
        assert(t_sxxi.range().size().x == _sizeX);
        assert(t_sxxi.range().size().y == _sizeY);
        assert(t_sxyi.range().size().x == _sizeX);
        assert(t_sxyi.range().size().y == _sizeY);
        assert(t_syyi.range().size().x == _sizeX);
        assert(t_syyi.range().size().y == _sizeY);
        assert(t_vxi.range().size().x == _sizeX);
        assert(t_vxi.range().size().y == _sizeY);
        assert(t_vyi.range().size().x == _sizeX);
        assert(t_vyi.range().size().y == _sizeY);
        assert(t_w1i.range().size().x == _sizeX);
        assert(t_w1i.range().size().y == _sizeY);
        assert(t_w1i_dy_p1_.range().size().x == _sizeX);
        assert(t_w1i_dy_p1_.range().size().y == _sizeY);
        assert(t_w1i_dy_p2_.range().size().x == _sizeX);
        assert(t_w1i_dy_p2_.range().size().y == _sizeY);
        assert(t_w1i_dy_m1_.range().size().x == _sizeX);
        assert(t_w1i_dy_m1_.range().size().y == _sizeY);
        assert(t_w1i_dy_m2_.range().size().x == _sizeX);
        assert(t_w1i_dy_m2_.range().size().y == _sizeY);
        assert(t_w2i.range().size().x == _sizeX);
        assert(t_w2i.range().size().y == _sizeY);
        assert(t_w2i_dy_p1_.range().size().x == _sizeX);
        assert(t_w2i_dy_p1_.range().size().y == _sizeY);
        assert(t_w2i_dy_p2_.range().size().x == _sizeX);
        assert(t_w2i_dy_p2_.range().size().y == _sizeY);
        assert(t_w2i_dy_m1_.range().size().x == _sizeX);
        assert(t_w2i_dy_m1_.range().size().y == _sizeY);
        assert(t_w2i_dy_m2_.range().size().x == _sizeX);
        assert(t_w2i_dy_m2_.range().size().y == _sizeY);
        assert(t_w3i.range().size().x == _sizeX);
        assert(t_w3i.range().size().y == _sizeY);
        assert(t_w3i_dy_p1_.range().size().x == _sizeX);
        assert(t_w3i_dy_p1_.range().size().y == _sizeY);
        assert(t_w3i_dy_p2_.range().size().x == _sizeX);
        assert(t_w3i_dy_p2_.range().size().y == _sizeY);
        assert(t_w3i_dy_m1_.range().size().x == _sizeX);
        assert(t_w3i_dy_m1_.range().size().y == _sizeY);
        assert(t_w3i_dy_m2_.range().size().x == _sizeX);
        assert(t_w3i_dy_m2_.range().size().y == _sizeY);
        assert(t_w4i.range().size().x == _sizeX);
        assert(t_w4i.range().size().y == _sizeY);
        assert(t_w4i_dy_p1_.range().size().x == _sizeX);
        assert(t_w4i_dy_p1_.range().size().y == _sizeY);
        assert(t_w4i_dy_p2_.range().size().x == _sizeX);
        assert(t_w4i_dy_p2_.range().size().y == _sizeY);
        assert(t_w4i_dy_m1_.range().size().x == _sizeX);
        assert(t_w4i_dy_m1_.range().size().y == _sizeY);
        assert(t_w4i_dy_m2_.range().size().x == _sizeX);
        assert(t_w4i_dy_m2_.range().size().y == _sizeY);
        auto _evaluateY2 = [&](Int _i, Int _j, auto type)
            __attribute__((always_inline,flatten))
        {
        using T = decltype(type);
        // temporaries
        T _c1;
        T _c2;
        T _dw1;
        T _dw2;
        T _dw3;
        T _dw4;
        // variables
        T _cpi;
        T _csi;
        T _lai;
        T _lmi;
        T _mui;
        T _sxxi;
        T _sxyi;
        T _syyi;
        T _vxi;
        T _vyi;
        T _w1i;
        T _w1i_dy_p1_;
        T _w1i_dy_p2_;
        T _w1i_dy_m1_;
        T _w1i_dy_m2_;
        T _w2i;
        T _w2i_dy_p1_;
        T _w2i_dy_p2_;
        T _w2i_dy_m1_;
        T _w2i_dy_m2_;
        T _w3i;
        T _w3i_dy_p1_;
        T _w3i_dy_p2_;
        T _w3i_dy_m1_;
        T _w3i_dy_m2_;
        T _w4i;
        T _w4i_dy_p1_;
        T _w4i_dy_p2_;
        T _w4i_dy_m1_;
        T _w4i_dy_m2_;
        // loads
        loadFromPtr(_cpi, &t_cpi.val(_i, _j));
        loadFromPtr(_csi, &t_csi.val(_i, _j));
        loadFromPtr(_lai, &t_lai.val(_i, _j));
        loadFromPtr(_lmi, &t_lmi.val(_i, _j));
        loadFromPtr(_mui, &t_mui.val(_i, _j));
        loadFromPtr(_sxxi, &t_sxxi.val(_i, _j));
        loadFromPtr(_sxyi, &t_sxyi.val(_i, _j));
        loadFromPtr(_syyi, &t_syyi.val(_i, _j));
        loadFromPtr(_vxi, &t_vxi.val(_i, _j));
        loadFromPtr(_vyi, &t_vyi.val(_i, _j));
        loadFromPtr(_w1i, &t_w1i.val(_i, _j));
        loadFromPtr(_w1i_dy_p1_, &t_w1i_dy_p1_.val(_i, _j));
        loadFromPtr(_w1i_dy_p2_, &t_w1i_dy_p2_.val(_i, _j));
        loadFromPtr(_w1i_dy_m1_, &t_w1i_dy_m1_.val(_i, _j));
        loadFromPtr(_w1i_dy_m2_, &t_w1i_dy_m2_.val(_i, _j));
        loadFromPtr(_w2i, &t_w2i.val(_i, _j));
        loadFromPtr(_w2i_dy_p1_, &t_w2i_dy_p1_.val(_i, _j));
        loadFromPtr(_w2i_dy_p2_, &t_w2i_dy_p2_.val(_i, _j));
        loadFromPtr(_w2i_dy_m1_, &t_w2i_dy_m1_.val(_i, _j));
        loadFromPtr(_w2i_dy_m2_, &t_w2i_dy_m2_.val(_i, _j));
        loadFromPtr(_w3i, &t_w3i.val(_i, _j));
        loadFromPtr(_w3i_dy_p1_, &t_w3i_dy_p1_.val(_i, _j));
        loadFromPtr(_w3i_dy_p2_, &t_w3i_dy_p2_.val(_i, _j));
        loadFromPtr(_w3i_dy_m1_, &t_w3i_dy_m1_.val(_i, _j));
        loadFromPtr(_w3i_dy_m2_, &t_w3i_dy_m2_.val(_i, _j));
        loadFromPtr(_w4i, &t_w4i.val(_i, _j));
        loadFromPtr(_w4i_dy_p1_, &t_w4i_dy_p1_.val(_i, _j));
        loadFromPtr(_w4i_dy_p2_, &t_w4i_dy_p2_.val(_i, _j));
        loadFromPtr(_w4i_dy_m1_, &t_w4i_dy_m1_.val(_i, _j));
        loadFromPtr(_w4i_dy_m2_, &t_w4i_dy_m2_.val(_i, _j));
        _c1 = _cpi * dt / dy;
        _c2 = _csi * dt / dy;

        _dw1 = advection(_c1, _w1i_dy_m2_, _w1i_dy_m1_, _w1i, _w1i_dy_p1_, _w1i_dy_p2_);
        _dw2 = advection(_c1, _w2i_dy_p2_, _w2i_dy_p1_, _w2i, _w2i_dy_m1_, _w2i_dy_m2_);
        _dw3 = advection(_c2, _w3i_dy_m2_, _w3i_dy_m1_, _w3i, _w3i_dy_p1_, _w3i_dy_p2_);
        _dw4 = advection(_c2, _w4i_dy_p2_, _w4i_dy_p1_, _w4i, _w4i_dy_m1_, _w4i_dy_m2_);

        // from_omega_y
        _vxi += _csi * (_dw3 - _dw4);
        _vyi += _cpi * (_dw1 - _dw2);
        _sxxi += _lai * (_dw1 + _dw2);
        _syyi += _lmi * (_dw1 + _dw2);
        _sxyi += _mui * (_dw3 + _dw4);
        storeToPtr(&t_sxxi.val(_i, _j),_sxxi);
        storeToPtr(&t_sxyi.val(_i, _j),_sxyi);
        storeToPtr(&t_syyi.val(_i, _j),_syyi);
        storeToPtr(&t_vxi.val(_i, _j),_vxi);
        storeToPtr(&t_vyi.val(_i, _j),_vyi);
        };
//        for (Int _j = 0; _j < _sizeY; _j++) {
//        Int _i = 0;
//        for (; _i < _sizeX / 4 * 4; _i+=4) {
//        _evaluateY2(_i, _j, Real4());
//        }
//        for (; _i < _sizeX; _i++) {
//        _evaluateY2(_i, _j, Real());
//        }
//        }

        struct Idx {
            Int x{}, y{}, z{};
        };
        Int xx = (_sizeX / 4);
        Int yy = (_sizeY / 2);
        Int tt = (time_steps * 4);
//        std::vector<Idx> indices;
//        indices.reserve(xx * yy * tt);

        auto timeStart = std::chrono::system_clock::now();

        // 33 seconds
//        tiledLoops(xx, yy, tt, [&] (int i, int j, int t) {
//            indices.push_back(Idx{i, j, t});
//        });

        // 31.5 seconds (half tile size: 24)
//        fixedTiledLoops(xx, yy, tt, 24, [&] (int i, int j, int t) {
//            indices.push_back(Idx{i, j, t});
//        });

        // 42 seconds
//        for (int t = 0; t < tt; t++)
//            for (int j = 0; j < yy; j++)
//                for (int i = 0; i < xx; i++)
//                    indices.push_back(Idx{i, j, t});

        auto timeEnd = std::chrono::system_clock::now();

        std::chrono::duration<double> timeDiff = timeEnd - timeStart;

        std::cerr << "Index gen time: " << timeDiff.count() << std::endl;

        timeStart = std::chrono::system_clock::now();


//        for (const Idx& idx : indices) {
//            Int i = idx.x;
//            Int j = idx.y;
//            Int t = idx.z;
//            switch(t % 4) {
//            case 0:
//                _evaluateX1(4*i, 2*j, Real4());
//                _evaluateX1(4*i, 2*j+1, Real4());
//                break;
//            case 1:
//                _evaluateX2(4*i, 2*j, Real4());
//                _evaluateX2(4*i, 2*j+1, Real4());
//                break;
//            case 2:
//                _evaluateY1(4*i, 2*j, Real4());
//                _evaluateY1(4*i, 2*j+1, Real4());
//                break;
//            case 3:
//                _evaluateY2(4*i, 2*j, Real4());
//                _evaluateY2(4*i, 2*j+1, Real4());
//                break;
//            }
//        }

        // 31.8 seconds
//        fixedTiledLoops(xx, yy, tt, 24,
//                   [&_evaluateX1,&_evaluateX2,&_evaluateY1,&_evaluateY2]
//                   (int i, int j, int t) __attribute__((always_inline)) {
//            switch(t % 4) {
//            case 0:
//                _evaluateX1(4*i, 2*j, Real4());
//                _evaluateX1(4*i, 2*j+1, Real4());
//                break;
//            case 1:
//                _evaluateX2(4*i, 2*j, Real4());
//                _evaluateX2(4*i, 2*j+1, Real4());
//                break;
//            case 2:
//                _evaluateY1(4*i, 2*j, Real4());
//                _evaluateY1(4*i, 2*j+1, Real4());
//                break;
//            case 3:
//                _evaluateY2(4*i, 2*j, Real4());
//                _evaluateY2(4*i, 2*j+1, Real4());
//                break;
//            }
//        });

        // 28.6 seconds (28.6775 nb)
//        tiledLoops(xx, yy, tt,
//                   [&_evaluateX1,&_evaluateX2,&_evaluateY1,&_evaluateY2]
//                   (int i, int j, int t) __attribute__((noinline)) {
//            switch(t % 4) {
//            case 0:
//                _evaluateX1(4*i, 2*j, Real4());
//                _evaluateX1(4*i, 2*j+1, Real4());
//                break;
//            case 1:
//                _evaluateX2(4*i, 2*j, Real4());
//                _evaluateX2(4*i, 2*j+1, Real4());
//                break;
//            case 2:
//                _evaluateY1(4*i, 2*j, Real4());
//                _evaluateY1(4*i, 2*j+1, Real4());
//                break;
//            case 3:
//                _evaluateY2(4*i, 2*j, Real4());
//                _evaluateY2(4*i, 2*j+1, Real4());
//                break;
//            }
//        });

        // 29.7 seconds
        //classicLoops(xx, yy, tt,
        //           [&_evaluateX1,&_evaluateX2,&_evaluateY1,&_evaluateY2]
        //           (int i, int j, int t) __attribute__((always_inline)) {
        //    switch(t % 4) {
        //    case 0:
        //        _evaluateX1(4*i, 2*j, Real4());
        //        _evaluateX1(4*i, 2*j+1, Real4());
        //        break;
        //    case 1:
        //        _evaluateX2(4*i, 2*j, Real4());
        //        _evaluateX2(4*i, 2*j+1, Real4());
        //        break;
        //    case 2:
        //        _evaluateY1(4*i, 2*j, Real4());
        //        _evaluateY1(4*i, 2*j+1, Real4());
        //        break;
        //    case 3:
        //        _evaluateY2(4*i, 2*j, Real4());
        //        _evaluateY2(4*i, 2*j+1, Real4());
        //        break;
        //    }
        //});

        // 43.86 nb
        tiledLoopsTemplated(xx, yy, tt,
                            [&_evaluateX1,&_evaluateX2,&_evaluateY1,&_evaluateY2]
                            (int i, int j, int t)
                            __attribute__((always_inline))
                            //__attribute__((noinline))
        {
            switch(t % 4) {
            case 0:
                _evaluateX1(4*i, 2*j, Real4());
                _evaluateX1(4*i, 2*j+1, Real4());
                break;
            case 1:
                _evaluateX2(4*i, 2*j, Real4());
                _evaluateX2(4*i, 2*j+1, Real4());
                break;
            case 2:
                _evaluateY1(4*i, 2*j, Real4());
                _evaluateY1(4*i, 2*j+1, Real4());
                break;
            case 3:
                _evaluateY2(4*i, 2*j, Real4());
                _evaluateY2(4*i, 2*j+1, Real4());
                break;
            }
        });

        // 29 seconds
//        for (Int t = 0; t < time_steps * 4; t++) {
//            for (Int j = 0; j < _sizeY; j++) {
//                Int i = 0;
//                for (; i < _sizeX / 4 * 4; i+=4) {
//                    switch(t % 4) {
//                    case 0:
//                        _evaluateX1(i, j, Real4());
//                        break;
//                    case 1:
//                        _evaluateX2(i, j, Real4());
//                        break;
//                    case 2:
//                        _evaluateY1(i, j, Real4());
//                        break;
//                    case 3:
//                        _evaluateY2(i, j, Real4());
//                        break;
//                    }
//                }
//                for (; i < _sizeX; i++) {
//                    switch(t % 4) {
//                    case 0:
//                        _evaluateX1(i, j, Real());
//                        break;
//                    case 1:
//                        _evaluateX2(i, j, Real());
//                        break;
//                    case 2:
//                        _evaluateY1(i, j, Real());
//                        break;
//                    case 3:
//                        _evaluateY2(i, j, Real());
//                        break;
//                    }
//                }
//            }
//        }

        timeEnd = std::chrono::system_clock::now();

        timeDiff = timeEnd - timeStart;

        std::cerr << "Looping time: " << timeDiff.count() << std::endl;

    }

    vxa.save(save_path + "/vx_" + withLeadingZeros(time_steps, 6) + ".bin");
    vya.save(save_path + "/vy_" + withLeadingZeros(time_steps, 6) + ".bin");
    sxxa.save(save_path + "/sxx_" + withLeadingZeros(time_steps, 6) + ".bin");
    syya.save(save_path + "/syy_" + withLeadingZeros(time_steps, 6) + ".bin");
    sxya.save(save_path + "/sxy_" + withLeadingZeros(time_steps, 6) + ".bin");

    std::cout << std::endl;

    return 0;
}
