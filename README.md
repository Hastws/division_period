# division_period

计算有理数 $M / N$ 十进制展开中的**循环节**，输出格式为 `integer.non_repeat[repeat]`。

项目包含两个可执行程序：

| 程序 | 说明 |
|------|------|
| `division_period` | 完整数论版本（欧拉函数求阶 + 多线程小数位计算） |
| `simple` | 简洁长除法版本（余数追踪法，适合教学理解） |

---

## 1. 数学原理

设约分后分数为：

$$\frac{M}{N}, \quad N = 2^{s_1} \cdot 5^{s_2} \cdot t, \quad \gcd(t, 10) = 1$$

| 量 | 含义 | 公式 |
|----|------|------|
| $p$ | 非循环部分长度 | $p = \max(s_1, s_2)$ |
| $\text{ord}_t(10)$ | 循环节长度 | 满足 $10^k \equiv 1 \pmod{t}$ 的最小正整数 $k$ |
| $q$ | 循环结束位置 | $q = p + \text{ord}_t(10)$ |

### 求阶算法

1. 对 $t$ 做质因数分解，计算欧拉函数 $\varphi(t)$
2. 分解 $\varphi(t)$ 的质因子
3. 由 Lagrange 定理，$\text{ord}_t(10) \mid \varphi(t)$
4. 从 $\varphi(t)$ 出发，逐个质因子尝试缩减指数，每次验证 $10^e \bmod t = 1$
5. 最终得到的 $e$ 即为 $\text{ord}_t(10)$

---

## 2. 代码结构

```
.
├── division_period.hpp   # 核心实现（header-only）
│   ├── Pow / PowMod      #   快速幂 / 模幂
│   ├── GCD               #   最大公因数（迭代）
│   ├── PrimePow          #   质因子幂 p^k
│   ├── Factor            #   质因数分解（含乘除/回滚）
│   ├── RemoveFactors2And5#   提取因子 2、5
│   ├── EulerPhi          #   欧拉函数
│   ├── CalOrd10T         #   求 ord_t(10)
│   └── DivisionPeriod    #   主计算类
├── normal.cpp            # 数论版主程序（内置测试 + CLI）
├── simple.cpp            # 长除法版主程序（CLI）
├── CMakeLists.txt        # 构建脚本
├── LICENSE
└── README.md
```

---

## 3. 构建与运行

### 3.1 环境要求

- CMake >= 3.16
- 支持 C++20 的编译器（GCC 10+ / Clang 12+）

### 3.2 构建

```bash
cmake -S . -B build
cmake --build build -j
```

### 3.3 运行

```bash
# 数论版 —— 运行内置测试用例
./build/division_period

# 数论版 —— 命令行指定分子分母
./build/division_period 1 7
./build/division_period 29 14
./build/division_period 1 7 1 896 100 5    # 一次计算多组

# 简洁版 —— 默认 1/6
./build/simple

# 简洁版 —— 命令行
./build/simple 22 7
```

---

## 4. 输出格式与示例

格式：`M / N == integer.non_repeat[repeat]`

```text
===== Built-in test cases =====
1 / 7 == 0.[142857]
1 / 896 == 0.0011160[714285]
1 / 1792 == 0.00055803[571428]
1 / 3584 == 0.000279017[857142]
1 / 7168 == 0.0001395089[285714]
1 / 14336 == 0.00006975446[428571]
1 / 1 == 1.[0]
11 / 1 == 11.[0]
100 / 5 == 20.[0]
8 / 512 == 0.015625[0]
4 / 512 == 0.0078125[0]
2 / 512 == 0.00390625[0]
1 / 512 == 0.001953125[0]
29 / 14 == 2.0[714285]
0 / 4 == 0.[0]
1 / 1000000007 == 0.[0000000009999999930000000489999996570000024009999831930001176489991764570057648009596463932824752470...]
```

- `[]` 内是循环节
- 有限小数的循环节记为 `[0]`
- 循环节超过 100 位时以 `...` 截断

---

## 5. 复杂度与实现细节

| 操作 | 复杂度 |
|------|--------|
| 质因数分解 | $O(\sqrt{n})$ |
| 模幂计算 | $O(\log k)$ |
| 小数位生成 | $O(q / 8)$，每 chunk 含 8 位十进制 |

当 chunk 数 >= 100 时，`DivisionPeriod` 自动将区间均分给 16 个线程并行计算。

---

## 6. 更新日志

### v2.0（当前版本）

1. **Bug Fix**: `Pow(base, 0)` 返回值修正为 `1`（原实现返回 `base`）
2. **Bug Fix**: 消除 `std::set` 元素通过 `const_cast` 修改导致的**未定义行为**
3. **优化**: 质因数分解循环改为 `k*k <= N`，避免重复 `sqrt` 调用
4. **优化**: `GCD` 改为迭代版，避免深递归栈溢出
5. **重构**: `GetNnm()` 重命名为 `GetNum()`，`Factorization()` 重命名为 `Factorize()`
6. **重构**: `Mod2and5()` 重命名为 `RemoveFactors2And5()`
7. **重构**: `#define MAX_OUTPUT` 改为 `constexpr`；错误输出走 `std::cerr`
8. **功能**: 两个程序均支持命令行参数输入
9. **功能**: `simple.cpp` 完全重写，逻辑清晰可读
10. **文档**: 全部源码添加完整 Doxygen 风格注释
11. **构建**: CMake 改用 `Threads::Threads`，添加 `-Wall -Wextra -Wpedantic`

---

## 7. 注意事项

- 分母必须为正整数（不可为 0）
- 内部使用 `uint32_t`，超范围乘法会被检测并报错
- `simple` 版本为单线程逐位计算，大分母下较慢

---

## 8. 许可证

本项目采用仓库中的 [LICENSE](LICENSE) 许可证。
