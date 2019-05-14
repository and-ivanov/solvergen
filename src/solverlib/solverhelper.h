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

template <typename T> struct isSequence { static const bool value = false; };

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

    Real& val(Int x, Int y) {
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

    Real& val(Int x, Int y) {
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
