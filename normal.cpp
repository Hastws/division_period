/**
 * @file   normal.cpp
 * @brief  division_period 主程序 —— 数论法计算有理数循环小数
 *
 * 用法：
 *   ./division_period              运行内置测试用例
 *   ./division_period M N          计算 M / N 的循环小数
 *   ./division_period M1 N1 M2 N2  一次计算多组
 *
 * 内置测试用例覆盖了：
 *   - 纯循环小数    (1/7)
 *   - 混合循环小数  (1/896, 29/14 ...)
 *   - 整除/有限小数 (1/1, 100/5, 8/512 ...)
 *   - 分子为 0      (0/4)
 *   - 大素数分母    (1/1000000007)
 */

#include "division_period.hpp"

#include <cstdlib>   // std::strtoul

int main(int argc, char *argv[]) {
    // ---- 命令行模式：M N [M N ...] ----
    if (argc >= 3 && (argc % 2 == 1)) {
        for (int i = 1; i + 1 < argc; i += 2) {
            const auto M = static_cast<uint32_t>(std::strtoul(argv[i],     nullptr, 10));
            const auto N = static_cast<uint32_t>(std::strtoul(argv[i + 1], nullptr, 10));
            if (N == 0) {
                std::cerr << "error: denominator must not be 0" << std::endl;
                return 1;
            }
            DivisionPeriod dp(M, N);
            dp.DispAllCyclicDecimal();
        }
        return 0;
    }

    if (argc != 1) {
        std::cerr << "usage: " << argv[0] << " [M N [M N ...]]" << std::endl;
        return 1;
    }

    // ---- 默认模式：运行内置测试用例 ----
    const std::vector<DivisionPeriod> cases = {
        // --- 纯循环小数：分母不含 2 和 5 的因子 ---
        DivisionPeriod(1, 7),       // 0.[142857]

        // --- 混合循环小数：分母含 2^s 或 5^s ---
        DivisionPeriod(1, 896),     // 0.0011160[714285]      896 = 2^7 * 7
        DivisionPeriod(1, 1792),    // 0.00055803[571428]     1792 = 2^8 * 7
        DivisionPeriod(1, 3584),    // 0.000279017[857142]    3584 = 2^9 * 7
        DivisionPeriod(1, 7168),    // 0.0001395089[285714]   7168 = 2^10 * 7
        DivisionPeriod(1, 14336),   // 0.00006975446[428571]  14336 = 2^11 * 7

        // --- 整除 / 有限小数（循环节为 [0]）---
        DivisionPeriod(1, 1),       // 1.[0]
        DivisionPeriod(11, 1),      // 11.[0]
        DivisionPeriod(100, 5),     // 20.[0]
        DivisionPeriod(8, 512),     // 0.015625[0]            512 = 2^9
        DivisionPeriod(4, 512),     // 0.0078125[0]
        DivisionPeriod(2, 512),     // 0.00390625[0]
        DivisionPeriod(1, 512),     // 0.001953125[0]

        // --- 分子大于分母 ---
        DivisionPeriod(29, 14),     // 2.0[714285]

        // --- 分子为 0 ---
        DivisionPeriod(0, 4),       // 0.[0]

        // --- 大素数分母（循环节很长，输出被截断）---
        DivisionPeriod(1, 1000000007),
        // 0.[0000000009999999930000000489999996570000024009999831930001176489991764570057648009596463932824752470...]
    };

    std::cout << "===== Built-in test cases =====" << std::endl;
    for (const auto &dp : cases) {
        dp.DispAllCyclicDecimal();
    }

    return 0;
}
