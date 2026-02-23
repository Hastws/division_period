/**
 * @file   simple.cpp
 * @brief  简洁版长除法 —— 用余数追踪法检测循环节
 *
 * 这是一个教学性的实现，展示了最基本的循环小数检测思路：
 *   1. 每步做长除法：remainder *= 10，取商位，取余
 *   2. 用 map 记录每个余数第一次出现的位置
 *   3. 一旦余数重复出现，两次之间的数字即为循环节
 *
 * 用法：
 *   ./simple              默认计算 1 / 6
 *   ./simple M N          计算 M / N
 */

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <vector>

int main(int argc, char *argv[]) {
    // ---- 解析参数 ----
    uint32_t numerator   = 1;   // 默认分子
    uint32_t denominator = 6;   // 默认分母

    if (argc == 3) {
        numerator   = static_cast<uint32_t>(std::strtoul(argv[1], nullptr, 10));
        denominator = static_cast<uint32_t>(std::strtoul(argv[2], nullptr, 10));
    } else if (argc != 1) {
        std::cerr << "usage: " << argv[0] << " [M N]" << std::endl;
        return 1;
    }

    if (denominator == 0) {
        std::cerr << "error: denominator must not be 0" << std::endl;
        return 1;
    }

    // ---- 长除法 ----
    std::vector<uint8_t> digits;                              // 小数各位
    std::unordered_map<uint32_t, std::size_t> first_seen;     // 余数 -> 首次出现位置

    uint32_t remainder = numerator % denominator;
    std::size_t repeat_start = 0;   // 循环节开始位置

    while (true) {
        // 余数为 0 => 整除，无循环节
        if (remainder == 0) {
            repeat_start = digits.size();
            break;
        }

        // 余数重复出现 => 从 first_seen 到当前位置为循环节
        if (const auto it = first_seen.find(remainder); it != first_seen.end()) {
            repeat_start = it->second;
            break;
        }

        // 记录该余数首次出现的位置，然后做一步长除
        first_seen[remainder] = digits.size();
        remainder *= 10;
        digits.push_back(static_cast<uint8_t>(remainder / denominator));
        remainder %= denominator;
    }

    // ---- 输出 ----
    std::cout << numerator << " / " << denominator << " = "
              << numerator / denominator << ".";

    // 非循环部分
    for (std::size_t i = 0; i < repeat_start; ++i) {
        std::cout << static_cast<int>(digits[i]);
    }

    // 循环部分
    std::cout << '[';
    if (repeat_start < digits.size()) {
        for (std::size_t i = repeat_start; i < digits.size(); ++i) {
            std::cout << static_cast<int>(digits[i]);
        }
    } else {
        std::cout << '0';   // 整除时显示 [0]
    }
    std::cout << ']' << std::endl;

    return 0;
}
