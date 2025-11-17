#ifndef DIVISION_PERIOD_AUTOALG_DIVISION_PERIOD_HPP
#define DIVISION_PERIOD_AUTOALG_DIVISION_PERIOD_HPP

#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>
#include <thread>
#include <set>
#include <cmath>
#include <limits>
#include <tuple>

#define MAX_OUTPUT  100

constexpr uint64_t MAX_INT32 = std::numeric_limits<uint32_t>::max();

// base^pow
inline uint32_t Pow(const uint32_t base, const uint8_t pow) {
    uint32_t res = base;
    for (uint8_t i = 1; i < pow; ++i) {
        res *= base;
    }
    return res;
}

// base^pow mod mod
inline uint32_t PowMod(uint64_t base, uint32_t pow, const uint32_t mod) {
    uint64_t res = 1;
    while (pow > 0) {
        if (pow % 2 == 1) {
            res = res * base % mod;
        }
        base = base * base % mod;
        pow = pow / 2;
    }
    return static_cast<uint32_t>(res);
}

// prime^pow
class PrimePow {
public:
    PrimePow(const PrimePow &obj) = default;

    PrimePow(const uint32_t prime, const uint8_t pow)
        : prime_(prime), pow_(pow) {
    }

    [[nodiscard]] uint32_t ToInt() const { return Pow(prime_, pow_); }
    [[nodiscard]] uint32_t GetPrime() const { return prime_; }
    [[nodiscard]] uint8_t GetPow() const { return pow_; }

    void SetPrime(const uint32_t pow) { prime_ = pow; }

protected:
    void TryMultiple(const PrimePow &pp) {
        pow_ += pp.pow_;
    }

    bool TryDiv(const PrimePow &pp) {
        if (prime_ != pp.prime_ || pow_ < pp.pow_) {
            return false;
        }

        pow_ -= pp.pow_;
        return (pow_ != 0);
    }

    friend bool operator <(const PrimePow &lhs, const PrimePow &rhs) {
        return lhs.prime_ < rhs.prime_;
    }

    friend std::ostream &operator <<(std::ostream &os, const PrimePow &pp) {
        os << "(" << pp.prime_ << "^" << static_cast<int>(pp.pow_) << ")";
        return os;
    }

    friend class Factor;

private:
    uint32_t prime_;
    uint8_t pow_;
};

// 质因数分解
class Factor {
public:
    Factor() = default;

    Factor(const Factor &obj) = default;

    Factor(const uint32_t num) {
        assert(num != 0);
        num_ = num;
        Factorization();
    }

    void DispCalResult() const {
        std::cout << num_ << " == " << *this << std::endl;
    }

    [[nodiscard]] uint32_t GetNnm() const { return num_; }
    [[nodiscard]] const std::set<PrimePow> &GetPrimePow() const { return factor_; }

    void SetNum(const uint32_t num) {
        assert(num != 0);
        num_ = num;
        factor_.clear();
        Factorization();
    }

    void Factorization() {
        if (!factor_.empty()) {
            factor_.clear();
        }

        uint32_t N = num_;
        auto sqrtN = static_cast<uint32_t>(std::sqrt(N));
        uint8_t counter = 0;

        // 单独处理因子 2
        while (N % 2 == 0) {
            N /= 2;
            ++counter;
        }
        if (counter != 0) {
            factor_.emplace(2, counter);
            sqrtN = static_cast<uint32_t>(std::sqrt(N));
        }

        // 枚举奇数因子
        counter = 0;
        for (uint32_t k = 3; k <= sqrtN; k += 2) {
            if (N % k == 0) {
                while (N % k == 0) {
                    N /= k;
                    ++counter;
                }
                factor_.emplace(k, counter);
                sqrtN = static_cast<uint32_t>(std::sqrt(N));
                counter = 0;
            }
        }

        // 若剩下的 N > 1，则本身是一个素数
        if (N != 1) {
            factor_.emplace(N, 1);
        }
    }

    // 尝试乘以一个 PrimePow
    bool TryMultiple(const PrimePow &prime_pow) {
        if (static_cast<uint64_t>(num_) *
            static_cast<uint64_t>(prime_pow.ToInt()) > MAX_INT32) {
            std::cout << "TryMultiple overflow, N=" << num_ << ", pp=" << prime_pow.ToInt() << std::endl;
            return false;
        }

        num_ *= prime_pow.ToInt();

        if (auto [it, is_success] = factor_.insert(prime_pow); !is_success) {
            auto &lhs = const_cast<PrimePow &>(*it);
            lhs.TryMultiple(prime_pow);
        }

        return true;
    }

    // 尝试除以一个 PrimePow
    bool TryDiv(const PrimePow &pp) {
        const auto it = factor_.find(pp);

        if (it == factor_.end() || it->pow_ < pp.pow_) {
            std::cout << "TryDiv dismatch, N=" << num_ << ", pp=" << pp.ToInt() << std::endl;
            return false;
        }

        num_ /= pp.ToInt();

        if (auto &lhs = const_cast<PrimePow &>(*it); !lhs.TryDiv(pp)) {
            factor_.erase(it);
        }

        return true;
    }

    bool TryMultiple(const Factor &ft) {
        for (auto it = ft.factor_.begin(); it != ft.factor_.end(); ++it) {
            if (!TryMultiple(*it)) {
                // 溢出，回滚之前的乘法
                for (auto it2 = ft.factor_.begin(); it2 != it; ++it2) { TryDiv(*it2); }
                std::cout << "TryMultiple overflow，N=" << num_ << ", ft.m_N=" << ft.num_ << std::endl;
                return false;
            }
        }
        return true;
    }

    // 尝试除以另一个 Factor
    bool TryDiv(const Factor &ft) {
        for (auto it = ft.factor_.begin(); it != ft.factor_.end(); ++it) {
            if (!TryDiv(*it)) {
                // 不匹配，回滚之前的除法
                for (auto it2 = ft.factor_.begin(); it2 != it; ++it2) {
                    TryMultiple(*it2);
                }
                std::cout << "TryDiv dismatch, N=" << num_
                        << ", ft.m_N=" << ft.num_ << std::endl;
                return false;
            }
        }
        return true;
    }

    // 输出操作符：打印为 (p1^s1)(p2^s2)...
    friend std::ostream &operator <<(std::ostream &os, const Factor &ft) {
        for (auto &pp: ft.factor_) {
            os << pp;
        }
        return os;
    }

private:
    uint32_t num_ = 1; // 原始整数 N
    std::set<PrimePow> factor_; // 质因数分解结果
};

/**
 * @brief 设 N = (2^s1)(5^s2)t
 *        除掉所有因子 2^s1 和 5^s2，返回引用 p = max{s1, s2}
 *        返回 Factor t
 */
inline Factor Mod2and5(const uint32_t N, uint32_t &p) {
    const Factor tmp_n(N); // 遍历用，禁止改变值，否则 for range 会炸
    Factor tmp_t(N); // 返回值

    uint8_t counter2 = 0;
    uint8_t counter5 = 0;
    uint32_t tmp_prime;

    // 除掉所有 2 的因子
    for (auto &pp: tmp_n.GetPrimePow()) {
        tmp_prime = pp.GetPrime();
        if (tmp_prime == 2) {
            counter2 = pp.GetPow();
            tmp_t.TryDiv(pp);
        }
    }

    // 除掉所有 5 的因子
    for (auto &pp: tmp_n.GetPrimePow()) {
        tmp_prime = pp.GetPrime();
        if (tmp_prime == 5) {
            counter5 = pp.GetPow();
            tmp_t.TryDiv(pp);
        }
    }

    p = static_cast<uint32_t>(counter2 > counter5 ? counter2 : counter5);
    return tmp_t;
}

/**
 * @brief 计算整数 T 的 EulerPhi(T)，返回一个 Factor 类型
 *        设 T = (p1^{s1}) * ... * (pn^{sn})
 *        EulerPhi(T) = (p1^{s1-1}) * ... * (pn^{sn-1}) * (p1-1) * ... * (pn-1)
 */
inline Factor EulerPhi(const Factor &t) {
    Factor tmpT(t);
    Factor tmpMult; // 空 Factor
    std::set<PrimePow> set_prime_pow; // 存 (p1^1) 到 (pn^1)

    // 把每个质因子取出来一个，放到 setPP
    for (auto pp: tmpT.GetPrimePow())
        set_prime_pow.emplace(pp.GetPrime(), 1);

    // 除掉每个 p_i
    for (auto &pp: set_prime_pow)
        tmpT.TryDiv(pp);

    // 乘以每个 (p_i - 1)
    for (auto &pp: set_prime_pow) {
        tmpMult.SetNum(pp.ToInt() - 1);
        tmpT.TryMultiple(tmpMult);
    }

    return tmpT;
}

/**
 * @brief 计算 ord_t 10
 *        设 EulerPhi(t) = (p1^s1)(p2^s2)...(pn^sn) = s
 *        对每个质因子 p_i： 尝试不断减少它的指数，检查 10^s (mod t) 是否仍为 1
 *        最后得到的指数组合即为 ord_t 10
 */
inline Factor CalOrd10T(const Factor &ft, const uint32_t mod) {
    Factor tmp_ft(ft);
    PrimePow tmpPP(1, 1); // 工具变量，m_prime 会变

    for (auto &pp: ft.GetPrimePow()) {
        tmpPP.SetPrime(pp.GetPrime()); // tmpPP = PrimePow(pp.GetPrime(), 1);
        uint8_t nowPow = pp.GetPow();

        while (nowPow--) {
            // 如果 nowPow 自减前为 0，说明该因子已经算完了
            tmp_ft.TryDiv(tmpPP);
            if (PowMod(10, tmp_ft.GetNnm(), mod) == 1) {
                // 继续尝试减少这个质因子的指数
                continue;
            }
            // 再减就不行了，恢复一次，换下一个质因子
            tmp_ft.TryMultiple(tmpPP);
            break;
        }
    }

    return tmp_ft;
}

/**
 * @brief 辗转相除法求最大公因子，输入可以有 0
 */
inline uint32_t GCD(const uint32_t a, const uint32_t b) {
    return b > 0 ? GCD(b, a % b) : a;
}

/**
 * @brief CycleDiv 类是计算循环除法的类，要求分母不为 0，输入参数非负
 *
 *        该类首先通过数论求阶的方法计算出 m_loopP 和 m_loopQ
 *        然后通过乘 10^n 法求余的办法求出每一位的十进制值
 */
class DivisionPeriod {
public:
    DivisionPeriod(uint32_t M, uint32_t N) : m_(M), n_(N) {
        assert(N != 0); // 分母不能为 0

        // 只需要计算小数部分（整数部分 m_M / m_N 直接输就好）
        M = M % N;

        const uint32_t tmpGCD = (M == 0 ? N : GCD(N, M));
        M /= tmpGCD; // 约分
        N /= tmpGCD;

        // 正式开始计算
        const Factor tmpT = Mod2and5(N, loop_p_); // N == (2^s1)(5^s2)(t)，m_loopP = max{s1,s2}
        const Factor tmpEP = EulerPhi(tmpT); // tmpEP = EulerPhi(t)
        loop_q_ = loop_p_ + CalOrd10T(tmpEP, tmpT.GetNnm()).GetNnm(); // 循环节终点

        const uint32_t loop_len = (loop_q_ / 8) + 1; // 每个元素存 8 位十进制


        decimal_.resize(loop_len); // 分配空间


        // 定义计算小数的 lambda 函数，方便多线程
        auto CalDecimal = [&](uint64_t MM, const uint64_t NN, uint32_t p, uint32_t q) {
            MM = MM * static_cast<uint64_t>(PowMod(100000000, p, NN)) % NN; // M = M * (10^(8p) % N)

            for (; p != q; ++p) {
                MM *= 100000000;
                decimal_[p] = static_cast<uint32_t>(MM / NN);
                MM %= NN;
            }
        };

        if (loop_len < 100) {
            // 比较小的就没必要用多线程加速了
            CalDecimal(M, N, 0, loop_len);
        } else {
            constexpr uint16_t thread_num = 16;

            std::vector<std::thread> threads;
            std::vector<uint32_t> p_and_q; // 装每个线程负责的 [p_i, q_i)

            for (uint16_t i = 0; i <= thread_num; ++i)
                p_and_q.emplace_back(
                    static_cast<uint32_t>(
                        static_cast<uint64_t>(loop_len) * static_cast<uint64_t>(i) / static_cast<uint64_t>(
                            thread_num)));

            threads.reserve(thread_num);
            for (uint16_t i = 0; i < thread_num; ++i)
                threads.emplace_back(CalDecimal, M, N, p_and_q[i], p_and_q[i + 1]);

            for (auto &tr: threads)
                tr.join();
        }
    }

    void DispAllCyclicDecimal() const {
        std::cout << m_ << " / " << n_ << " == " << *this;
    }

    [[nodiscard]] uint32_t GetM() const { return m_; }
    [[nodiscard]] uint32_t GetN() const { return n_; }
    [[nodiscard]] uint32_t GetP() const { return loop_p_; } // 非循环部分长度
    [[nodiscard]] uint32_t GetQ() const { return loop_q_; } // 循环结束位置
    [[nodiscard]] uint32_t GetLoopLength() const { return loop_q_ - loop_p_; } // 循环节长度
    [[nodiscard]] uint32_t GetDecimal(const uint32_t N) const { return decimal_[N]; }

    // 输出操作符
    friend std::ostream &operator <<(std::ostream &os, const DivisionPeriod &cd);

private:
    uint32_t m_; // 分子（原始）
    uint32_t n_; // 分母（原始）

    std::vector<uint32_t> decimal_; // 小数部分，每个储存 8 位十进制小数
    uint32_t loop_p_{}; // 循环节范围 (m_loopP, m_loopQ]
    uint32_t loop_q_;
};

inline std::ostream &operator <<(std::ostream &os, const DivisionPeriod &cd) {
    uint32_t tmp_p = cd.loop_p_;
    uint32_t tmp_q = cd.loop_q_ > MAX_OUTPUT ? MAX_OUTPUT : cd.loop_q_; // 最大只显示 MAX_OUTPUT 位

    os << cd.m_ / cd.n_ << "."; // 整数部分

    uint32_t k = 0;

    // 非循环部分（前半段整 8 位）
    for (k = 0; k < tmp_p / 8; ++k)
        os << std::setw(8) << std::setfill('0') << static_cast<uint32_t>(cd.decimal_[k]);

    // 非循环部分（如果不是 8 的整倍数，补尾巴）
    if (tmp_p % 8)
        os << std::setw(tmp_p % 8) << std::setfill('0')
                << (cd.decimal_[k] / Pow(10, 8 - tmp_p % 8));

    os << "["; // 循环部分开始

    if (tmp_p / 8 == tmp_q / 8) {
        // 循环部分在同一个 uint32 里
        os << std::setw(tmp_q - tmp_p) << std::setfill('0')
                << (cd.decimal_[k] / Pow(10, 8 - tmp_q % 8)
                    % Pow(10, tmp_q - tmp_p));
    } else {
        // 先输出当前这个 uint32 里属于循环的后半段
        os << std::setw(8 - tmp_p % 8) << std::setfill('0')
                << (cd.decimal_[k] % Pow(10, 8 - tmp_p % 8));

        // 中间完整的若干个 uint32
        for (++k; k < tmp_q / 8; ++k)
            os << std::setw(8) << std::setfill('0') << static_cast<uint32_t>(cd.decimal_[k]);

        // 最后一个不满 8 位的部分
        if (tmp_q % 8)
            os << std::setw(tmp_q % 8) << std::setfill('0')
                    << (cd.decimal_[k] / Pow(10, 8 - tmp_q % 8));
    }

    if (tmp_q != cd.loop_q_) // 如果被截断为了 MAX_OUTPUT，输出 ...
        os << "...";

    os << "]" << std::endl;
    return os;
}


#endif // DIVISION_PERIOD_AUTOALG_DIVISION_PERIOD_HPP
