import re
import sys


ATOMIC_WEIGHTS = {
    'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811,
    'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797,
    'Na': 22.98976928, 'Mg': 24.3050, 'Al': 26.9815386, 'Si': 28.0855, 'P': 30.973762,
    'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
    'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
    'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.64, 'As': 74.92160, 'Se': 78.96, 'Br': 79.904,
    'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585, 'Zr': 91.224,
    'Nb': 92.90638, 'Mo': 95.96, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.90550,
    'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.710,
    'Sb': 121.760, 'Te': 127.60, 'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519,
    'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116, 'Pr': 140.90765, 'Nd': 144.242,
    'Pm': 145.0, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535,
    'Dy': 162.500, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.054,
    'Lu': 174.9668, 'Hf': 178.49, 'Ta': 180.94788, 'W': 183.84, 'Re': 186.207,
    'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.592,
    'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.98040, 'Po': 209.0, 'At': 210.0,
    'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0, 'Ac': 227.0, 'Th': 232.03806,
    'Pa': 231.03588, 'U': 238.02891
}

BRACKET_PAIRS = {'(': ')', '[': ']'}
CLOSING_BRACKETS = set(BRACKET_PAIRS.values())


def _tokenize(formula):
    tokens = []
    i = 0
    while i < len(formula):
        c = formula[i]
        if c in BRACKET_PAIRS:
            tokens.append(('OPEN', c))
            i += 1
        elif c in CLOSING_BRACKETS:
            tokens.append(('CLOSE', c))
            i += 1
        elif c.isdigit():
            j = i
            while j < len(formula) and formula[j].isdigit():
                j += 1
            tokens.append(('NUM', int(formula[i:j])))
            i = j
        elif c.isupper():
            j = i + 1
            while j < len(formula) and formula[j].islower():
                j += 1
            elem = formula[i:j]
            tokens.append(('ELEM', elem))
            i = j
        else:
            raise ValueError(f"无法识别的字符: '{c}' (位置 {i})")
    return tokens


def _parse_tokens(tokens, pos):
    elements = {}
    while pos < len(tokens):
        kind, value = tokens[pos]
        if kind == 'OPEN':
            sub_elements, pos = _parse_tokens(tokens, pos + 1)
            multiplier = 1
            if pos < len(tokens) and tokens[pos][0] == 'NUM':
                multiplier = tokens[pos][1]
                pos += 1
            for elem, count in sub_elements.items():
                elements[elem] = elements.get(elem, 0) + count * multiplier
        elif kind == 'CLOSE':
            return elements, pos + 1
        elif kind == 'ELEM':
            elem = value
            count = 1
            if pos + 1 < len(tokens) and tokens[pos + 1][0] == 'NUM':
                count = tokens[pos + 1][1]
                pos += 2
            else:
                pos += 1
            elements[elem] = elements.get(elem, 0) + count
        elif kind == 'NUM':
            pos += 1
        else:
            pos += 1
    return elements, pos


def parse_formula(formula):
    formula = formula.strip()
    if not formula:
        raise ValueError("化学式不能为空")
    tokens = _tokenize(formula)
    elements, _ = _parse_tokens(tokens, 0)
    return sorted(elements.items())


def calculate_molecular_weight(formula):
    elements = parse_formula(formula)
    if not elements:
        raise ValueError("化学式解析结果为空")
    total = 0.0
    for elem, count in elements:
        if elem not in ATOMIC_WEIGHTS:
            raise ValueError(f"未知元素: {elem}")
        total += ATOMIC_WEIGHTS[elem] * count
    return round(total, 2)


def main():
    if len(sys.argv) > 1:
        formulas = sys.argv[1:]
    else:
        print("分子量计算器 (输入化学式，如 H2O、Fe2(SO4)3、[Cu(NH3)4]SO4)")
        print("输入 q 退出\n")
        while True:
            try:
                formula = input("请输入化学式: ").strip()
            except (EOFError, KeyboardInterrupt):
                print()
                break
            if formula.lower() == 'q':
                break
            if not formula:
                continue
            try:
                elements = parse_formula(formula)
                expanded = " + ".join(f"{e}×{c}" for e, c in elements)
                weight = calculate_molecular_weight(formula)
                print(f"  展开: {expanded}")
                print(f"  分子量: {weight:.2f} g/mol\n")
            except ValueError as e:
                print(f"  错误: {e}\n")
        return

    for formula in formulas:
        try:
            weight = calculate_molecular_weight(formula)
            print(f"{formula}: {weight:.2f} g/mol")
        except ValueError as e:
            print(f"{formula}: 错误 - {e}")


if __name__ == '__main__':
    main()
