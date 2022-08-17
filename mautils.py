import math
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product, groupby
from math import gcd, lcm
from sympy.ntheory import factorint


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


def get_lcm_for_fractions(*fractions: Fraction, ignore_zero=True):
    if ignore_zero:
        fractions = [fraction for fraction in fractions if fraction != 0]
    fractions_gcd = gcd(*[fraction.denominator for fraction in fractions])
    fractions_lcm = Fraction(lcm(*[fraction.numerator for fraction in fractions]), fractions_gcd)
    return fractions_lcm


def get_lcm_for_combinations(fractions: set[Fraction, Fraction, ...]) -> list[Fraction]:
    lcm_combinations = set()
    for i in range(2, len(fractions) + 1):
        for combination in combinations(fractions, r=i):
            lcm_combinations.add(get_lcm_for_fractions(*combination))
    return sorted(lcm_combinations)


def get_possible_lcm_configurations_for_fractions(*fractions: Fraction):
    fractions = sorted(set(fractions))
    lcm_configurations = []
    lcm_values = []

    for reference_fraction in fractions:
        lcm_configurations.append([fraction / reference_fraction for fraction in fractions])
        lcm_values.append(get_lcm_for_fractions(*lcm_configurations[-1]))

    return tuple(fractions), tuple(lcm_values), tuple(lcm_configurations)


def get_inverted_fractions(*fractions: Fraction):
    return [Fraction(fraction.denominator, fraction.numerator) for fraction in fractions]


def get_all_12_tet_chords():
    # Chord defined as any set of octave reduced notes in root position which does not contain a semitone intervall
    all_key_combinations = product((0, 1), repeat=10)
    relevant_key_combinations = []
    for comb in all_key_combinations:
        for i in range(len(comb) - 1):
            if comb[i] == comb[i + 1] == 1:
                break
        else:
            relevant_key_combinations.append((1, 0, *comb))

    return relevant_key_combinations


def chord_to_fractions(chord):
    if len(chord) != 12:
        raise ValueError(f"chord must be len 12, was {len(chord)}")
    keyboard_as_fractions = [Fraction(1, 1), Fraction(16, 15), Fraction(9, 8), Fraction(6, 5), Fraction(5, 4),
                             Fraction(4, 3), Fraction(7, 5), Fraction(3, 2), Fraction(8, 5), Fraction(5, 3),
                             Fraction(16, 9), Fraction(15, 8)]
    fractions = [0] * 12
    for i, key in enumerate(chord):
        if key:
            fractions[i] = keyboard_as_fractions[i]

    return fractions


def sort_fractions_by_lcm(fractions: list[list, ...]) -> dict[list[set[Fraction, ...]]]:
    lcm_sorted = defaultdict(list)
    for fraction_list in fractions:
        fraction_set = {fraction for fraction in fraction_list if fraction != 0}
        lcm_sorted[get_lcm_for_fractions(*fraction_set)].append(fraction_set)

    return dict(sorted(lcm_sorted.items()))


def get_fractions_with_bases_up_to(fractions: set[Fraction, ...], max_base=5):
    filtered_fractions = set()
    for fraction in fractions:
        for factor in factorint(fraction.numerator):
            if factor > max_base:
                break
        else:
            for factor in factorint(fraction.denominator):
                if factor > max_base:
                    break
            else:
                filtered_fractions.add(fraction)

    return filtered_fractions


def to_octave_reduced_fractions(fractions, octave_size=2):
    reduced_fractions = set()
    for fraction in fractions:
        while fraction < octave_size:
            fraction *= octave_size
        fraction /= octave_size
        reduced_fractions.add(fraction)
    return reduced_fractions
