// Author is Marek Cerny (c) 2024, me@marekcerny.com

#pragma once 


#include <Eigen/Core>

#include <array>
#include <cmath>
#include <compare>
#include <iostream>
#include <string>




namespace wl {

namespace crt {
    using ull = unsigned long long;

    // prime source:) https://www.math.utah.edu/~pa/MDS/primes.html

    constexpr auto primes = std::array{
        2147483647ull,	
        2147483629ull,	
        2147483587ull,	
        2147483579ull, 
        2147483563ull,	
        2147483549ull,	
        2147483543ull, 
        2147483497ull, 
        2147483489ull,	
        2147483477ull,	
        2147483423ull, 
        2147483399ull,
        2147483353ull,
        2147483323ull, 
        2147483269ull, 
        2147483249ull,
        2147483237ull,
        2147483179ull,
        2147483171ull,
        2147483137ull,
    };

    template<ull R>
    struct crtu64t {
        static_assert(R < primes.size(), "R must be less than primes.size()");

        crtu64t() {
            for (auto& mod : mods) mod = 0ull;
        }

        crtu64t(ull val) {
            for (int r = 0; r < R; ++r) {
                mods[r] = val % primes[r];
            }
        }

        static auto min() -> crtu64t {
            return crtu64t(0ull);
        }

        static auto max() -> crtu64t {
            crtu64t res;
            for (int r = 0; r < R; ++r) {
                res.mods[r] = primes[r] - 1;
            }
            return res;
        }

        auto operator+=(crtu64t const& other) -> crtu64t& {
            for (int r = 0; r < R; ++r) {
                mods[r] = (mods[r] + other.mods[r]) % primes[r];
            }
            return *this;
        }

        friend auto operator+(crtu64t lhs, crtu64t const& rhs) -> crtu64t {
            lhs += rhs;
            return lhs;
        }

        auto operator*=(crtu64t const& other) -> crtu64t& {
            for (int r = 0; r < R; ++r) {
                mods[r] = (mods[r] * other.mods[r]) % primes[r];
            }
            return *this;
        }

        friend auto operator*(crtu64t lhs, crtu64t const& rhs) -> crtu64t {
            lhs *= rhs;
            return lhs;
        }

        auto operator-=(crtu64t const& other) -> crtu64t& {
            for (int r = 0; r < R; ++r) {
                mods[r] = (mods[r] + primes[r] - other.mods[r]) % primes[r];
            }
            return *this;
        }

        friend auto operator-(crtu64t lhs, crtu64t const& rhs) -> crtu64t {
            lhs -= rhs;
            return lhs;
        }

        auto operator-() const -> crtu64t {
            return crtu64t(0ULL) - *this;
        }

        auto operator<=>(crtu64t const& other) const -> std::strong_ordering {
            for (int r = 0; r < R; ++r) {
                if (mods[r] != other.mods[r]) {
                    return mods[r] <=> other.mods[r];
                }
            }
            return std::strong_ordering::equal;
        }

        auto operator==(crtu64t const& other) const -> bool {
            return (*this <=> other) == std::strong_ordering::equal;
        }

        friend auto to_string(crtu64t const& val) -> std::string {
            std::string res;
            for (int r = 0; r < R; ++r) {
                res += std::to_string(val.mods[r]);
                if (r != R-1) res += "x_";
            }
            return res;
        }

        friend auto operator<<(std::ostream& os, crtu64t const& val) -> std::ostream& {
            for (int r = 0; r < R; ++r) {
                os << val.mods[r];
                if (r != R-1) os << "x_";
            }
            return os;
        }

        friend auto log10(crtu64t const& val) -> double {
            double res = 0.0;
            for (int r = 0; r < R; ++r) {
                res += std::log10(static_cast<double>(val.mods[r]));
            }
            return res;
        }


        std::array<ull, R> mods;
    };


}

} // namespace wl


namespace Eigen {
    template<std::size_t R>
    struct NumTraits<wl::crt::crtu64t<R>> : GenericNumTraits<wl::crt::crtu64t<R>> {
        typedef wl::crt::crtu64t<R> Real;

        enum {
            IsComplex = 0,
            IsInteger = 1,
            IsSigned = 0,
            RequireInitialization = 1,
            ReadCost = 10 * R,
            AddCost = 10 * R,
            MulCost = 20 * R
        };

        // // Since it's an integral type, epsilon and precision are not applicable
        static inline Real epsilon() {
            return Real(0);
        }

        static inline Real dummy_precision() {
            return Real(0);
        }

        // For an integer type, digits10 can represent the maximum representable value
        static inline int digits10() {
            double log10_max = log10(Real::max());
            return static_cast<int>(std::ceil(log10_max));
        }

        static inline Real highest() {
            return Real::max();
        }

        static inline Real lowest() {
            return Real::min();
        }
    };
} // namespace Eigen