import re
import sys


ATOMIC_WEIGHTS = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.81,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
    'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904,
    'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
    'Nb': 92.906, 'Mo': 95.95, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.906,
    'Pd': 106.42, 'Ag': 107.868, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.710,
    'Sb': 121.760, 'Te': 127.60, 'I': 126.904, 'Xe': 131.293, 'Cs': 132.905,
    'Ba': 137.327, 'La': 138.905, 'Ce': 140.116, 'Pr': 140.908, 'Nd': 144.242,
    'Pm': 145.0, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.925,
    'Dy': 162.500, 'Ho': 164.930, 'Er': 167.259, 'Tm': 168.934, 'Yb': 173.045,
    'Lu': 174.967, 'Hf': 178.486, 'Ta': 180.948, 'W': 183.84, 'Re': 186.207,
    'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084, 'Au': 196.967, 'Hg': 200.592,
    'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.980, 'Po': 209.0, 'At': 210.0,
    'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0, 'Ac': 227.0, 'Th': 232.038,
    'Pa': 231.036, 'U': 238.029
}


def parse_formula(formula):
    pattern = r'([A-Z][a-z]?)(\d*)'
    tokens = []
    pos = 0
    while pos < len(formula):
        if formula[pos] == '(':
            depth = 1
            pos += 1
            start = pos
            while pos < len(formula) and depth > 0:
                if formula[pos] == '(':
                    depth += 1
                elif formula[pos] == ')':
                    depth -= 1
                pos += 1
            inner = formula[start:pos - 1]
            multiplier = ''
            while pos < len(formula) and formula[pos].isdigit():
                multiplier += formula[pos]
                pos += 1
            mult = int(multiplier) if multiplier else 1
            inner_tokens = parse_formula(inner)
            for elem, count in inner_tokens:
                tokens.append((elem, count * mult))
        else:
            match = re.match(pattern, formula[pos:])
            if match and match.group(1):
                elem = match.group(1)
                count_str = match.group(2)
                count = int(count_str) if count_str else 1
                tokens.append((elem, count))
                pos += len(match.group(0))
            else:
                raise ValueError(f"无法解析化学式中的字符: '{formula[pos]}' 在位置 {pos}")
    return tokens


def calculate_molecular_weight(formula):
    tokens = parse_formula(formula)
    total_weight = 0.0
    for elem, count in tokens:
        if elem not in ATOMIC_WEIGHTS:
            raise ValueError(f"未知元素: {elem}")
        total_weight += ATOMIC_WEIGHTS[elem] * count
    return round(total_weight, 2)


def main():
    if len(sys.argv) > 1:
        formulas = sys.argv[1:]
    else:
        formulas = ['H2O', 'C6H12O6', 'NaCl', 'Ca(OH)2', 'H2SO4', 'CO2']
        print("未输入化学式，使用默认示例:\n")

    for formula in formulas:
        try:
            weight = calculate_molecular_weight(formula)
            print(f"{formula}: {weight:.2f} g/mol")
        except ValueError as e:
            print(f"{formula}: 错误 - {e}")


if __name__ == '__main__':
    main()
