from math import factorial

TOL = 1e-10


def difference_table(y):
    table = [y[:]]

    while len(table[-1]) > 1:
        previous = table[-1]
        current = [
            previous[i + 1] - previous[i]
            for i in range(len(previous) - 1)
        ]
        table.append(current)

    return table


def print_table(x, table):
    print("\nDifference Table")

    print(f"{'x':>10}", end="")
    for k in range(len(table)):
        print(f"{'f(x)' if k == 0 else f'D{k}f':>12}", end="")
    print()

    for i in range(len(x)):
        print(f"{x[i]:>10g}", end="")

        for column in table:
            if i < len(column):
                print(f"{column[i]:>12g}", end="")
            else:
                print(f"{'':>12}", end="")
        print()


def stirling(table, m, p):
    result = table[0][m]

    for k in range(1, len(table)):
        if k % 2 == 0:
            r = k // 2
            coefficient = p * p

            for j in range(1, r):
                coefficient *= p * p - j * j

            row = m - r

            if 0 <= row < len(table[k]):
                result += coefficient * table[k][row] / factorial(k)

        else:
            r = (k - 1) // 2
            coefficient = p

            for j in range(1, r + 1):
                coefficient *= p * p - j * j

            row1 = m - r - 1
            row2 = m - r

            if 0 <= row1 < len(table[k]) and 0 <= row2 < len(table[k]):
                result += (
                    coefficient
                    * (table[k][row1] + table[k][row2])
                    / (2 * factorial(k))
                )

    return result


def bessel(table, m, p):
    m = m - 1
    result = table[0][m] + p * table[1][m]

    for k in range(2, len(table)):
        if k % 2 == 0:
            r = k // 2
            coefficient = p * (p - 1)

            for j in range(1, r):
                coefficient *= (p + j) * (p - j - 1)

            row = m - r

            if 0 <= row and row + 1 < len(table[k]):
                result += (
                    coefficient
                    * (table[k][row] + table[k][row + 1])
                    / (2 * factorial(k))
                )

        else:
            r = (k - 1) // 2
            coefficient = p * (p - 1) * (p - 0.5)

            for j in range(1, r):
                coefficient *= (p + j) * (p - j - 0.5)

            row = m - r

            if 0 <= row < len(table[k]):
                result += coefficient * table[k][row] / factorial(k)

    return result


def everett(table, m, p):
    m = m - 1
    q = 1 - p

    result = q * table[0][m] + p * table[0][m + 1]

    for k in range(2, len(table), 2):
        r = k // 2

        q_coefficient = q * (q * q - 1)
        p_coefficient = p * (p * p - 1)

        for j in range(2, r + 1):
            q_coefficient *= q * q - j * j
            p_coefficient *= p * p - j * j

        row = m - r

        if 0 <= row and row + 1 < len(table[k]):
            result += (
                q_coefficient * table[k][row]
                + p_coefficient * table[k][row + 1]
            ) / factorial(k + 1)

    return result


def gaussian_forward(table, m, p):
    result = table[0][m]

    for k in range(1, len(table)):
        if k % 2 == 1:
            r = (k - 1) // 2
            coefficient = p

            for j in range(1, r + 1):
                coefficient *= (p + j) * (p - j)

            row = m - r

        else:
            r = k // 2
            coefficient = p

            for j in range(1, r + 1):
                coefficient *= p - j

            for j in range(1, r):
                coefficient *= p + j

            row = m - r

        if 0 <= row < len(table[k]):
            result += coefficient * table[k][row] / factorial(k)

    return result


def gaussian_backward(table, m, p):
    result = table[0][m]

    for k in range(1, len(table)):
        if k % 2 == 1:
            r = (k - 1) // 2
            coefficient = p

            for j in range(1, r + 1):
                coefficient *= (p + j) * (p - j)

            row = m - r - 1

        else:
            r = k // 2
            coefficient = p

            for j in range(1, r + 1):
                coefficient *= p + j

            for j in range(1, r):
                coefficient *= p - j

            row = m - r

        if 0 <= row < len(table[k]):
            result += coefficient * table[k][row] / factorial(k)

    return result


print("=" * 70)
print("INTERPOLATION WITH CENTRAL DIFFERENCES")
print("=" * 70)

n = int(input("Enter number of data points : "))

if n < 3:
    print("At least 3 data points are required.")
    raise SystemExit

x = list(map(float, input("x values : ").split()))
y = list(map(float, input("f(x) values : ").split()))

if len(x) != n or len(y) != n:
    print("Error: Number of x and f(x) values must be equal to n.")
    raise SystemExit

h = x[1] - x[0]

if h <= 0:
    print("Error: x values must be in increasing order.")
    raise SystemExit

for i in range(1, n - 1):
    if abs((x[i + 1] - x[i]) - h) > TOL:
        print("Error: x values must be equally spaced.")
        raise SystemExit

table = difference_table(y)

print_table(x, table)

print("\nAvailable Formulae")
print("1. Stirling")
print("2. Bessel")
print("3. Everett")
print("4. Gaussian Forward")
print("5. Gaussian Backward")

choice = int(input("\nChoose a formula (1-5) : "))

if choice < 1 or choice > 5:
    print("Invalid choice.")
    raise SystemExit

query = float(input("Enter value of x for interpolation : "))

if query < x[0] - TOL or query > x[-1] + TOL:
    print("Query is outside the tabulated range.")
    raise SystemExit

for i in range(n):
    if abs(query - x[i]) < TOL:
        print(f"\nf({query:g}) = {y[i]:g}  (tabulated value)")
        raise SystemExit

m = n // 2

p = (query - x[m]) / h

results = {}

results["Stirling"] = stirling(table, m, p)

bessel_m = m - 1
bessel_p = (query - x[bessel_m]) / h
results["Bessel"] = bessel(table, m, bessel_p)

everett_p = bessel_p
results["Everett"] = everett(table, m, everett_p)

results["Gaussian Forward"] = gaussian_forward(table, m, p)
results["Gaussian Backward"] = gaussian_backward(table, m, p)

formula_names = {
    1: "Stirling",
    2: "Bessel",
    3: "Everett",
    4: "Gaussian Forward",
    5: "Gaussian Backward"
}

selected = formula_names[choice]

print("\nSelected Formula")
print("-" * 50)
print(f"Formula : {selected}")
print(f"x       : {query:g}")
print(f"Result  : {results[selected]:.10g}")

print("\nOverall Comparison")
print("-" * 50)

for name, result in results.items():
    print(f"{name:<25}: {result:.10g}")

values = list(results.values())
average = sum(values) / len(values)

print("-" * 50)
print(f"{'Average of all results':<25}: {average:.10g}")