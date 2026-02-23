/**
 * @file   division_period.hpp
 * @brief  计算有理数 M/N 十进制展开的循环节
 *
 * 核心数学原理：
 *   对于有理数 M/N，先约分得 M'/N'，再把 N' 分解为
 *       N' = 2^s1 * 5^s2 * t,  其中 gcd(t, 10) = 1
 *   非循环部分长度 p = max(s1, s2)
 *   循环节长度 = ord_t(10)  即 10 在模 t 意义下的乘法阶
 *   循环结束位置 q = p + ord_t(10)
 *
 * 求阶算法：
 *   1. 计算 EulerPhi(t) 并分解质因子
 *   2. 从 EulerPhi(t) 出发，逐个质因子尝试缩减指数
 *      每次检查 10^(缩减后的值) mod t 是否等于 1
 *   3. 最终得到的值就是 ord_t(10)
 *
 * 小数位计算：
 *   将小数每 8 位打包进一个 uint32_t（相当于 10^8 进制）
 *   对较大循环节可多线程并行计算各段
 */

#ifndef DIVISION_PERIOD_HPP
#define DIVISION_PERIOD_HPP

#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <thread>
#include <vector>

/// 输出时最多展示的小数位数（超出部分以 ... 截断）
constexpr uint32_t MAX_OUTPUT = 100;

/// uint32_t 的最大值，用于溢出检查
constexpr uint64_t MAX_UINT32 = std::numeric_limits<uint32_t>::max();

// ============================================================================
//  基础数学工具
// ============================================================================

/**
 * @brief 整数快速幂  base^exp
 * @note  不检查溢出，调用者需保证结果在 uint32_t 范围内
 */
inline uint32_t Pow(const uint32_t base, const uint8_t exp) {
    uint32_t res = 1;
    for (uint8_t i = 0; i < exp; ++i) {
        res *= base;
    }
    return res;
}

/**
 * @brief 模幂  base^exp mod mod，使用快速幂（二进制分解）
 * @note  中间结果用 uint64_t 防溢出
 */
inline uint32_t PowMod(uint64_t base, uint32_t exp, const uint32_t mod) {
    uint64_t res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) {                     // exp 为奇数
            res = res * base % mod;
        }
        base = base * base % mod;          // base = base^2 mod mod
        exp >>= 1;
    }
    return static_cast<uint32_t>(res);
}

/**
 * @brief 辗转相除法求最大公因数（迭代版，避免深递归栈溢出）
 */
inline uint32_t GCD(uint32_t a, uint32_t b) {
    while (b > 0) {
        const uint32_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

// ============================================================================
//  PrimePow —— 表示单个质因子的幂  p^k
// ============================================================================

/**
 * @brief 表示 prime^pow 这样一个质因子幂
 *
 * 按 prime 排序（< 比较），用于存入 std::set<PrimePow>
 * TryMultiple / TryDiv 仅供 Factor 类通过友元调用
 */
class PrimePow {
public:
    PrimePow(const PrimePow &) = default;
    PrimePow(const uint32_t prime, const uint8_t pow)
        : prime_(prime), pow_(pow) {}

    /// 还原为整数值 prime^pow
    [[nodiscard]] uint32_t ToInt()     const { return Pow(prime_, pow_); }
    [[nodiscard]] uint32_t GetPrime()  const { return prime_; }
    [[nodiscard]] uint8_t  GetPow()    const { return pow_; }

    /// 修改底数（CalOrd10T 中复用同一对象时需要）
    void SetPrime(const uint32_t prime) { prime_ = prime; }

protected:
    /// 指数相加：this->pow_ += pp.pow_
    void TryMultiple(const PrimePow &pp) {
        pow_ += pp.pow_;
    }

    /**
     * @brief 指数相减：this->pow_ -= pp.pow_
     * @return 减完后 pow_ 是否仍 > 0（若为 0 表示该因子应被移除）
     */
    bool TryDiv(const PrimePow &pp) {
        if (prime_ != pp.prime_ || pow_ < pp.pow_) {
            return false;
        }
        pow_ -= pp.pow_;
        return (pow_ != 0);
    }

    /// 按质数大小排序
    friend bool operator<(const PrimePow &lhs, const PrimePow &rhs) {
        return lhs.prime_ < rhs.prime_;
    }

    /// 格式化输出 (prime^pow)
    friend std::ostream &operator<<(std::ostream &os, const PrimePow &pp) {
        os << "(" << pp.prime_ << "^" << static_cast<int>(pp.pow_) << ")";
        return os;
    }

    friend class Factor;

private:
    uint32_t prime_;   ///< 质数底
    uint8_t  pow_;     ///< 指数
};

// ============================================================================
//  Factor —— 整数的质因数分解
// ============================================================================

/**
 * @brief 将正整数分解为质因数之积，内部用 std::set<PrimePow> 有序存储
 *
 * 提供 TryMultiple / TryDiv 来动态修改分解结果（含溢出/不匹配检查与回滚）
 * 欧拉函数、求阶等算法都基于此类操作。
 */
class Factor {
public:
    Factor() = default;
    Factor(const Factor &) = default;

    /// 从正整数构造，自动做质因数分解
    explicit Factor(const uint32_t num) {
        assert(num != 0);
        num_ = num;
        Factorize();
    }

    /// 打印分解结果，如  12 == (2^2)(3^1)
    void DispCalResult() const {
        std::cout << num_ << " == " << *this << std::endl;
    }

    [[nodiscard]] uint32_t                   GetNum()      const { return num_; }
    [[nodiscard]] const std::set<PrimePow> & GetPrimePow() const { return factor_; }

    /// 重置为新的正整数并重新分解
    void SetNum(const uint32_t num) {
        assert(num != 0);
        num_ = num;
        factor_.clear();
        Factorize();
    }

    /**
     * @brief 质因数分解（试除法）
     *
     * 先提取所有因子 2，再枚举奇数 k，判断 k*k <= N
     * 最后若 N > 1 则 N 本身是一个质数
     */
    void Factorize() {
        factor_.clear();

        uint32_t N = num_;
        uint8_t counter = 0;

        // ---- 提取因子 2 ----
        while (N % 2 == 0) {
            N /= 2;
            ++counter;
        }
        if (counter != 0) {
            factor_.emplace(2, counter);
        }

        // ---- 枚举奇数因子 ----
        counter = 0;
        for (uint32_t k = 3;
             static_cast<uint64_t>(k) * k <= N;   // 用 uint64_t 防 k*k 溢出
             k += 2)
        {
            if (N % k == 0) {
                while (N % k == 0) {
                    N /= k;
                    ++counter;
                }
                factor_.emplace(k, counter);
                counter = 0;
            }
        }

        // ---- 剩下的 N 若 >1，则自身为质数 ----
        if (N != 1) {
            factor_.emplace(N, 1);
        }
    }

    // ----------------------------------------------------------------
    //  乘以 / 除以单个 PrimePow
    // ----------------------------------------------------------------

    /**
     * @brief 将当前因数乘以 prime_pow
     * @return 成功返回 true；若会导致 uint32_t 溢出则返回 false
     */
    bool TryMultiple(const PrimePow &prime_pow) {
        if (static_cast<uint64_t>(num_) * prime_pow.ToInt() > MAX_UINT32) {
            std::cerr << "[TryMultiple] overflow: N=" << num_
                      << ", pp=" << prime_pow.ToInt() << std::endl;
            return false;
        }

        num_ *= prime_pow.ToInt();

        // 如果 factor_ 中已有同底质因子，合并指数；否则直接插入
        const auto it = factor_.find(prime_pow);
        if (it == factor_.end()) {
            factor_.insert(prime_pow);
        } else {
            PrimePow merged = *it;
            merged.TryMultiple(prime_pow);
            factor_.erase(it);
            factor_.insert(merged);
        }
        return true;
    }

    /**
     * @brief 将当前因数除以 pp
     * @return 成功返回 true；若因子不匹配或指数不足则返回 false
     */
    bool TryDiv(const PrimePow &pp) {
        const auto it = factor_.find(pp);
        if (it == factor_.end() || it->pow_ < pp.pow_) {
            std::cerr << "[TryDiv] mismatch: N=" << num_
                      << ", pp=" << pp.ToInt() << std::endl;
            return false;
        }

        num_ /= pp.ToInt();

        // 拷贝-删除-重插入，避免修改 set 元素的 UB
        PrimePow reduced = *it;
        factor_.erase(it);
        if (reduced.TryDiv(pp)) {       // pow_ 仍 > 0，重新插入
            factor_.insert(reduced);
        }
        return true;
    }

    // ----------------------------------------------------------------
    //  乘以 / 除以另一个 Factor（逐质因子操作，失败时回滚）
    // ----------------------------------------------------------------

    /**
     * @brief 逐质因子乘以 ft；若某步溢出则回滚之前所有乘法
     */
    bool TryMultiple(const Factor &ft) {
        for (auto it = ft.factor_.begin(); it != ft.factor_.end(); ++it) {
            if (!TryMultiple(*it)) {
                // 回滚已完成的乘法
                for (auto it2 = ft.factor_.begin(); it2 != it; ++it2) {
                    TryDiv(*it2);
                }
                std::cerr << "[TryMultiple] overflow: N=" << num_
                          << ", ft=" << ft.num_ << std::endl;
                return false;
            }
        }
        return true;
    }

    /**
     * @brief 逐质因子除以 ft；若某步不匹配则回滚之前所有除法
     */
    bool TryDiv(const Factor &ft) {
        for (auto it = ft.factor_.begin(); it != ft.factor_.end(); ++it) {
            if (!TryDiv(*it)) {
                for (auto it2 = ft.factor_.begin(); it2 != it; ++it2) {
                    TryMultiple(*it2);
                }
                std::cerr << "[TryDiv] mismatch: N=" << num_
                          << ", ft=" << ft.num_ << std::endl;
                return false;
            }
        }
        return true;
    }

    /// 打印为 (p1^s1)(p2^s2)... 形式
    friend std::ostream &operator<<(std::ostream &os, const Factor &ft) {
        for (const auto &pp : ft.factor_) {
            os << pp;
        }
        return os;
    }

private:
    uint32_t num_ = 1;            ///< 当前整数值
    std::set<PrimePow> factor_;   ///< 质因数分解结果（按质数升序）
};

// ============================================================================
//  数论辅助函数
// ============================================================================

/**
 * @brief 从 N 中提取所有因子 2 和 5
 *
 * 设 N = 2^s1 * 5^s2 * t，其中 gcd(t, 10) == 1
 * @param[in]  N   待分解的正整数
 * @param[out] p   max(s1, s2)，即非循环部分长度
 * @return     Factor t（已去除 2 和 5 的因子）
 */
inline Factor RemoveFactors2And5(const uint32_t N, uint32_t &p) {
    const Factor original(N);   // 只读，用于遍历
    Factor t(N);                // 将从中除去 2、5

    uint8_t s1 = 0;   // 2 的指数
    uint8_t s2 = 0;   // 5 的指数

    for (const auto &pp : original.GetPrimePow()) {
        if (pp.GetPrime() == 2) {
            s1 = pp.GetPow();
            t.TryDiv(pp);
        } else if (pp.GetPrime() == 5) {
            s2 = pp.GetPow();
            t.TryDiv(pp);
        }
    }

    p = (s1 > s2) ? s1 : s2;
    return t;
}

/**
 * @brief 计算欧拉函数 EulerPhi(T) 并返回其质因数分解
 *
 * 公式：若 T = p1^a1 * p2^a2 * ... * pn^an
 *       则 phi(T) = p1^(a1-1)*(p1-1) * p2^(a2-1)*(p2-1) * ... * pn^(an-1)*(pn-1)
 *
 * 实现方式：
 *   1. 对每个质因子 pi，把 T 除以 pi（降一次指数）
 *   2. 再乘以 (pi - 1)
 */
inline Factor EulerPhi(const Factor &t) {
    Factor result(t);
    Factor temp;                          // 临时 Factor 用于乘法
    std::set<PrimePow> primes;            // 收集所有 (pi, 1)

    // 收集每个质因子
    for (const auto &pp : result.GetPrimePow()) {
        primes.emplace(pp.GetPrime(), 1);
    }

    // 除掉每个 pi（指数降 1）
    for (const auto &pp : primes) {
        result.TryDiv(pp);
    }

    // 乘以每个 (pi - 1)
    for (const auto &pp : primes) {
        temp.SetNum(pp.ToInt() - 1);      // temp = Factor(pi - 1)
        result.TryMultiple(temp);
    }

    return result;
}

/**
 * @brief 计算 ord_t(10) —— 10 在模 t 意义下的乘法阶
 *
 * 由 Lagrange 定理，ord_t(10) | phi(t)
 * 算法：从 phi(t) 的分解出发，逐个质因子尝试缩减指数
 *       每次检查 10^(当前值) mod t == 1?
 *       如果是则继续缩减，否则恢复并换下一个质因子
 *
 * @param phi_t   EulerPhi(t) 的质因数分解
 * @param t       模数（已去除 2 和 5 的因子）
 * @return        ord_t(10) 的质因数分解
 */
inline Factor CalOrd10T(const Factor &phi_t, const uint32_t t) {
    Factor ord(phi_t);
    PrimePow one_prime(1, 1);   // 复用对象，prime 会被替换

    for (const auto &pp : phi_t.GetPrimePow()) {
        one_prime.SetPrime(pp.GetPrime());
        uint8_t remaining = pp.GetPow();

        while (remaining--) {
            ord.TryDiv(one_prime);         // 尝试除掉一个 pi
            if (PowMod(10, ord.GetNum(), t) == 1) {
                continue;                  // 仍满足 10^ord == 1 (mod t)，继续缩减
            }
            ord.TryMultiple(one_prime);    // 不行了，恢复
            break;
        }
    }

    return ord;
}

// ============================================================================
//  DivisionPeriod —— 计算并展示有理数的循环小数
// ============================================================================

/**
 * @brief 计算 M/N 的十进制展开，并识别循环节
 *
 * 构造时即完成全部计算，随后可通过 DispAllCyclicDecimal() 或 operator<< 输出。
 * 小数部分以 8 位十进制为一组存入 vector<uint32_t>，
 * 当 chunk 数 >= 100 时自动启用 16 线程并行计算以加速。
 */
class DivisionPeriod {
public:
    /**
     * @brief 构造函数：完成所有计算
     * @param M 分子（非负）
     * @param N 分母（正整数）
     */
    DivisionPeriod(uint32_t M, uint32_t N) : m_(M), n_(N) {
        assert(N != 0);

        // 取小数部分的分子
        M = M % N;

        // 约分
        const uint32_t g = (M == 0) ? N : GCD(N, M);
        M /= g;
        N /= g;

        // ---- 核心数论计算 ----
        // 把 N 分解为 2^s1 * 5^s2 * t，得 p = max(s1, s2)
        const Factor t = RemoveFactors2And5(N, loop_p_);

        // ord_t(10) = 循环节长度
        const Factor phi_t = EulerPhi(t);
        const uint32_t ord  = CalOrd10T(phi_t, t.GetNum()).GetNum();
        loop_q_ = loop_p_ + ord;

        // ---- 小数位计算 ----
        // 每个 uint32_t 存 8 位十进制（即 10^8 进制）
        const uint32_t chunk_count = (loop_q_ / 8) + 1;
        decimal_.resize(chunk_count);

        // lambda：计算 [p, q) 范围内各 chunk 的值
        // MM = M * 10^(8*p) mod NN，然后逐步 *10^8 取商和余
        auto CalDecimal = [this](uint64_t MM, const uint64_t NN,
                                 uint32_t p, uint32_t q) {
            MM = MM * static_cast<uint64_t>(PowMod(100000000, p, static_cast<uint32_t>(NN))) % NN;
            for (; p != q; ++p) {
                MM *= 100000000ULL;
                decimal_[p] = static_cast<uint32_t>(MM / NN);
                MM %= NN;
            }
        };

        if (chunk_count < 100) {
            // 数据量小，单线程即可
            CalDecimal(M, N, 0, chunk_count);
        } else {
            // 数据量大，16 线程并行
            constexpr uint16_t kThreadNum = 16;

            // 计算各线程负责的 chunk 区间 [p_i, q_i)
            std::vector<uint32_t> bounds;
            bounds.reserve(kThreadNum + 1);
            for (uint16_t i = 0; i <= kThreadNum; ++i) {
                bounds.push_back(static_cast<uint32_t>(
                    static_cast<uint64_t>(chunk_count) * i / kThreadNum));
            }

            std::vector<std::thread> threads;
            threads.reserve(kThreadNum);
            for (uint16_t i = 0; i < kThreadNum; ++i) {
                threads.emplace_back(CalDecimal, M, N, bounds[i], bounds[i + 1]);
            }
            for (auto &thr : threads) {
                thr.join();
            }
        }
    }

    /// 打印完整的循环小数表示
    void DispAllCyclicDecimal() const {
        std::cout << m_ << " / " << n_ << " == " << *this;
    }

    // ---- 访问接口 ----
    [[nodiscard]] uint32_t GetM()          const { return m_; }
    [[nodiscard]] uint32_t GetN()          const { return n_; }
    [[nodiscard]] uint32_t GetP()          const { return loop_p_; }
    [[nodiscard]] uint32_t GetQ()          const { return loop_q_; }
    [[nodiscard]] uint32_t GetLoopLength() const { return loop_q_ - loop_p_; }
    [[nodiscard]] uint32_t GetDecimalChunk(const uint32_t idx) const {
        return decimal_[idx];
    }

    friend std::ostream &operator<<(std::ostream &os, const DivisionPeriod &cd);

private:
    uint32_t m_;   ///< 原始分子
    uint32_t n_;   ///< 原始分母

    std::vector<uint32_t> decimal_;   ///< 小数位（每元素 8 位十进制）
    uint32_t loop_p_{};               ///< 非循环部分长度
    uint32_t loop_q_{};               ///< 循环部分结束位置（即 p + 循环节长度）
};

/**
 * @brief 格式化输出 DivisionPeriod
 *
 * 格式：integer.non_repeat[repeat]
 * 超过 MAX_OUTPUT 位时在 ] 前加 ...
 */
inline std::ostream &operator<<(std::ostream &os, const DivisionPeriod &cd) {
    const uint32_t p = cd.loop_p_;
    const uint32_t q = (cd.loop_q_ > MAX_OUTPUT) ? MAX_OUTPUT : cd.loop_q_;

    // ---- 整数部分 ----
    os << cd.m_ / cd.n_ << ".";

    uint32_t k = 0;

    // ---- 非循环部分：完整的 8 位 chunk ----
    for (k = 0; k < p / 8; ++k) {
        os << std::setw(8) << std::setfill('0') << cd.decimal_[k];
    }

    // ---- 非循环部分：末尾不足 8 位的零头 ----
    if (p % 8) {
        os << std::setw(p % 8) << std::setfill('0')
           << (cd.decimal_[k] / Pow(10, 8 - p % 8));
    }

    // ---- 循环部分 ----
    os << "[";

    if (p / 8 == q / 8) {
        // 循环部分与非循环尾巴在同一个 chunk 内
        os << std::setw(q - p) << std::setfill('0')
           << (cd.decimal_[k] / Pow(10, 8 - q % 8) % Pow(10, q - p));
    } else {
        // 当前 chunk 属于循环的后半段
        os << std::setw(8 - p % 8) << std::setfill('0')
           << (cd.decimal_[k] % Pow(10, 8 - p % 8));

        // 中间完整的 chunk
        for (++k; k < q / 8; ++k) {
            os << std::setw(8) << std::setfill('0') << cd.decimal_[k];
        }

        // 最后一个 chunk 的前几位
        if (q % 8) {
            os << std::setw(q % 8) << std::setfill('0')
               << (cd.decimal_[k] / Pow(10, 8 - q % 8));
        }
    }

    if (q != cd.loop_q_) {
        os << "...";   // 被截断
    }

    os << "]" << std::endl;
    return os;
}

#endif // DIVISION_PERIOD_HPP
