#include "division_period.hpp"

using namespace std;

int main() {
    const std::vector vector = {
        DivisionPeriod(1, 7), // expected: 1 / 7 == 0.[142857]
        DivisionPeriod(1, 896), // expected: 1 / 896 == 0.0011160[714285]
        DivisionPeriod(1, 1792), // expected: 1 / 1792 == 0.00055803[571428]
        DivisionPeriod(1, 3584), // expected: 1 / 3584 == 0.000279017[857142]
        DivisionPeriod(1, 7168), // expected: 1 / 7168 == 0.0001395089[285714]
        DivisionPeriod(1, 14336), // expected: 1 / 14336 == 0.00006975446[428571]
        DivisionPeriod(1, 1), // expected: 1 / 1 == 1.[0]
        DivisionPeriod(29, 14), // expected: 29 / 14 == 2.0[714285]
        DivisionPeriod(0, 4), // expected: 0 / 14 == 0.[0]
        DivisionPeriod(11, 1), // expected: 11 / 1 == 11.[0]
        DivisionPeriod(100, 5), // expected: 100 / 5 == 20.[0]
        DivisionPeriod(8, 512), // expected: 8 / 512 == 0.015625[0]
        DivisionPeriod(4, 512), // expected: 4 / 512 == 0.0078125[0]
        DivisionPeriod(2, 512), // expected: 2 / 512 == 0.00390625[0]
        DivisionPeriod(1, 512), // expected: 1 / 512 == 0.001953125[0]
        DivisionPeriod(1, 1000000007)
        // expected: 1 / 1000000007 == 0.[0000000009999999930000000489999996570000024009999831930001176489991764570057648009596463932824752470...]
    };

    for (auto &division_period: vector) {
        division_period.DispAllCyclicDecimal();
    }
    return 0;
}
