#ifndef SOLVERHELPER_H
#define SOLVERHELPER_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <iomanip>
#include <cassert>

#ifndef NDEBUG
#include <fenv.h>
#endif

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

struct Size {
    explicit Size(Int x, Int y, Int z = 1) : x(x), y(y), z(z) {}
    Int x, y, z;
};

bool operator==(const Size& l, const Size& r) {
    return l.x == r.x && l.y == r.y && l.z == r.z;
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
    Range& z(Int s, Int e) {
        zs = s;
        ze = e;
        return *this;
    }
    Range dx(Int d) {
        return Range().x(xs + d, xe + d).y(ys, ye).z(zs, ze);
    }
    Range dy(Int d) {
        return Range().x(xs, xe).y(ys + d, ye + d).z(zs, ze);
    }
    Range dz(Int d) {
        return Range().x(xs, xe).y(ys, ye).z(zs + d, ze + d);
    }
    Size size() const { return Size(xe - xs, ye - ys, ze - zs); }

    Int xs = 0;
    Int xe = -1;
    Int ys = 0;
    Int ye = -1;
    Int zs = 0;
    Int ze = -1;
};



struct Array {
    explicit Array(Int nx, Int ny, Int nz=1)
        : Array(Size(nx, ny, nz)) {}
    explicit Array(Size size)
        : size(size)
        , data(size_t(size.x * size.y * size.z))
    {}
    Array& load(const std::string& s) {
        std::ifstream ifs(s);
        CHECK2(ifs.is_open(), "Can't open file for reading: " + s);
        ifs.seekg(0, std::ios_base::end);
        auto fileSize = ifs.tellg();
        ifs.seekg(0, std::ios_base::beg);
        Int arraySize = size.x * size.y * size.z;
        CHECK2(fileSize == arraySize * sizeof(Real),
               "Size of data in file '" + s + "' is not equal to requested" +
               "(" + std::to_string(size.x) +
               "," + std::to_string(size.y) +
               "," + std::to_string(size.z) + ")");
        ifs.read(reinterpret_cast<char*>(data.data()), arraySize * sizeof(Real));
        return *this;
    }
    Array& save(const std::string& s) {
        std::ofstream ofs(s);
        CHECK2(ofs.is_open(), "Can't open file for writing: " + s);
        Int arraySize = size.x * size.y * size.z;
        ofs.write(reinterpret_cast<char*>(data.data()), arraySize * sizeof(Real));
        return *this;
    }

    Real& val(Int x, Int y, Int z = 0) {
        return data[size_t(x + size.x * y + size.x * size.y * z)];
    }

    Size size;
    std::vector<Real> data;
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

template <typename T> struct isSequence { static const bool value = false; };

struct Expr;

struct Assign {};

struct Sequence {};

struct Var {
    explicit Var(Array& a, Range r = Range()) : array(a), range(r) {
        if (range.xe == -1) range.xe = a.size.x;
        if (range.ye == -1) range.ye = a.size.y;
        if (range.ze == -1) range.ze = a.size.z;
        assert(0 <= range.xs);
        assert(range.xs <= range.xe);
        assert(range.xe <= array.size.x);
        assert(0 <= range.ys);
        assert(range.ys <= range.ye);
        assert(range.ye <= array.size.y);
        assert(0 <= range.zs);
        assert(range.zs <= range.ze);
        assert(range.ze <= array.size.z);
    }
    Var& x(Int s, Int e) { range.xs = s; range.xe = e; return *this; }
    Var& y(Int s, Int e) { range.ys = s; range.ye = e; return *this; }
    Var& z(Int s, Int e) { range.zs = s; range.ze = e; return *this; }

    Var dx(Int d) { return Var(array, range.dx(d)); }
    Var dy(Int d) { return Var(array, range.dy(d)); }
    Var dz(Int d) { return Var(array, range.dz(d)); }

    Assign operator=(const Var&) { return Assign(); }
    template <typename Expr, std::enable_if_t<isExprNode<Expr>::value, int> = 0>
    Assign operator=(const Expr&) { return Assign(); }

    Assign operator+=(const Expr&) { return Assign(); }
    Assign operator-=(const Expr&) { return Assign(); }
    Assign operator*=(const Expr&) { return Assign(); }
    Assign operator/=(const Expr&) { return Assign(); }

    Real& val(Int x, Int y, Int z = 0) {
        return array.val(range.xs + x, range.ys + y, range.zs + z);
    }

    Array& array;
    Range range;
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
          std::enable_if_t<isExprNode<Expr1>::value, int> = 0,
          std::enable_if_t<isExprNode<Expr2>::value, int> = 0>
Expr operator+(Expr1, Expr2) { return Expr(); }

template <typename Expr1, typename Expr2,
          std::enable_if_t<isExprNode<Expr1>::value, int> = 0,
          std::enable_if_t<isExprNode<Expr2>::value, int> = 0>
Expr operator-(Expr1, Expr2) { return Expr(); }

template <typename Expr1, typename Expr2,
          std::enable_if_t<isExprNode<Expr1>::value, int> = 0,
          std::enable_if_t<isExprNode<Expr2>::value, int> = 0>
Expr operator*(Expr1, Expr2) { return Expr(); }

template <typename Expr1, typename Expr2,
          std::enable_if_t<isExprNode<Expr1>::value, int> = 0,
          std::enable_if_t<isExprNode<Expr2>::value, int> = 0>
Expr operator/(Expr1, Expr2) { return Expr(); }

template <typename E,
          std::enable_if_t<isExprNode<E>::value, int> = 0>
Expr operator+(E) { return Expr(); }

template <typename E,
          std::enable_if_t<isExprNode<E>::value, int> = 0>
Expr operator-(E) { return Expr(); }

template <typename E, typename X,
          std::enable_if_t<isExprNode<E>::value, int> = 0,
          std::enable_if_t<isExprNode<X>::value, int> = 0>
Expr operator+=(E, X) { return Expr(); }


template <typename Expr1, typename Expr2,
          std::enable_if_t<isExprNode<Expr1>::value, int> = 0,
          std::enable_if_t<isExprNode<Expr2>::value, int> = 0>
bool operator>(Expr1, Expr2) { return true; }

template <typename Expr1, typename Expr2,
          std::enable_if_t<isExprNode<Expr1>::value, int> = 0,
          std::enable_if_t<isExprNode<Expr2>::value, int> = 0>
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


#endif // SOLVERHELPER_H
