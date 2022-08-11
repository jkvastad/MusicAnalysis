import math
from fractions import Fraction
from itertools import combinations
from math import gcd, lcm


def get_closest_scientific_pitch(fraction):
    c_4 = 261.6256
    twelve_tet_names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    fraction_as_pitch = fraction * c_4
    smallest_diff = 100000
    name_number = '-'
    for n in range(-48, 60):
        scientific_pitch = c_4 * 2 ** (n / 12)
        current_diff = abs(fraction_as_pitch - scientific_pitch)
        if current_diff < smallest_diff:
            smallest_diff = current_diff
            name_number = n
    pitch_name = twelve_tet_names[name_number % 12] + f'{4 + name_number // 12}'
    closest_pitch = c_4 * 2 ** (name_number / 12)
    return fraction, pitch_name, 1200 * math.log2(
        fraction_as_pitch / closest_pitch), smallest_diff / closest_pitch, smallest_diff


def get_lcm_for_fractions(*fractions: Fraction):
    fractions_gcd = gcd(*[fraction.denominator for fraction in fractions])
    fractions_lcm = Fraction(lcm(*[fraction.numerator for fraction in fractions]), fractions_gcd)
    return fractions_lcm


def get_lcm_for_combinations(fractions: set[Fraction, Fraction, ...]) -> list[Fraction]:
    lcm_combinations = set()
    for i in range(2, len(fractions) + 1):
        for combination in combinations(fractions, r=i):
            lcm_combinations.add(get_lcm_for_fractions(*combination))
    return sorted(lcm_combinations)