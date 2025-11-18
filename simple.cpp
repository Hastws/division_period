#include <cstdint>
#include <iostream>

static int remainder = 1;
static constexpr int N = 6;

int main() {
    std::uint8_t decimal[N];
    int power[N];

    int len = 0;
    int p = 0;
    int q = 0;

    for (int &i: power) {
        i = N;
    }

    while (true) {
        decimal[len] = static_cast<std::uint8_t>(remainder / N);
        remainder = remainder % N;
        if (remainder == 0) {
            p = q = len;
            break;
        }
        if (power[remainder] != N) {
            p = power[remainder];
            q = len;
            break;
        }

        power[remainder] = len;
        ++len;
        remainder *= 10;
    }

    std::cout << remainder << " / " << N << " = " << 1 / 7 << ".";

    for (int i = 1; i <= p; ++i) {
        std::cout << static_cast<int>(decimal[i]);
    }
    std::cout << '[';
    if (p != q) {
        for (int i = p + 1; i <= q; ++i) {
            std::cout << static_cast<int>(decimal[i]);
        }
    }
    std::cout << ']';
    return 0;
}
