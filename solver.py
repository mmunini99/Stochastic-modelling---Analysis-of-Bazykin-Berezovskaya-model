import sympy as sp

# Definire le variabili simboliche
A, B, C, D = sp.symbols('A B C D')

# Impostare i parametri
r = 1.0
k = 1.0
l = 0.5
m = 0.78

# Definire le equazioni
eq1 = 2 * r * m * (-2 * m + l + k) * A - m * (B + C) + (r * m**2 - r * m * k)**2
eq2 = -r * m**2 * (2 * B + A) + r * l * m * (B + A) + r * k * m * (B + A) - m * D - r * l * k * A
eq3 = -r * m**2 * (2 * C + A) + r * l * m * (C + A) + r * k * m * (C + A) - m * D - r * l * k * A
eq4 = B + C - r * (m - l) * (k - m)

# Risolvere il sistema
solutions = sp.solve([eq1, eq2, eq3, eq4], (A, B, C, D))

# Visualizzare le soluzioni
print(solutions)