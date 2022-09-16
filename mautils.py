import math
from collections import defaultdict, deque
from fractions import Fraction
from itertools import combinations, product
from math import gcd, lcm
from sympy.ntheory import factorint
from enum import Enum


class Scale(Enum):
    MAJOR = [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1]
    HARMONIC_MAJOR = [1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1]
    HARMONIC_MINOR = [1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1]
    MELODIC_MINOR = [1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    OCTACTONIC = [1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0]
    BLUES = [1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0]


MAJOR_SCALE_FRACTIONS = {Fraction(1), Fraction(9, 8), Fraction(5, 4), Fraction(4, 3), Fraction(3, 2), Fraction(5, 3),
                         Fraction(15, 8)}
NOTE_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']


def get_scale_name(scale_to_identify: list) -> (list, str):
    scale_to_identify = deque(scale_to_identify)
    fundamental = 0
    for i in range(len(scale_to_identify)):
        for scale in Scale:
            if list(scale_to_identify) == scale.value:
                return scale.name, NOTE_NAMES[fundamental]
        scale_to_identify.rotate()
        fundamental -= 1


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


def get_lcm_for_fractions(fractions: set[Fraction, ...], ignore_zero=True) -> Fraction:
    if ignore_zero:
        fractions = [fraction for fraction in fractions if fraction != 0]
    fractions_gcd = gcd(*[fraction.denominator for fraction in fractions])
    fractions_lcm = Fraction(lcm(*[fraction.numerator for fraction in fractions]), fractions_gcd)
    return fractions_lcm


def get_lcm_for_combinations(fractions: set[Fraction, Fraction, ...]) -> list[Fraction]:
    lcm_combinations = set()
    for i in range(2, len(fractions) + 1):
        for combination in combinations(fractions, r=i):
            lcm_combinations.add(get_lcm_for_fractions(set(combination)))
    return sorted(lcm_combinations)


def get_possible_lcm_configurations_for_fractions(*fractions: Fraction):
    fractions = sorted(set(fractions))
    lcm_configurations = []
    lcm_values = []

    for reference_fraction in fractions:
        lcm_configurations.append([fraction / reference_fraction for fraction in fractions])
        lcm_values.append(get_lcm_for_fractions(*lcm_configurations[-1]))

    return tuple(fractions), tuple(lcm_values), tuple(lcm_configurations)


def get_inverted_fractions(*fractions: Fraction) -> list[Fraction, ...]:
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


def get_octave_reduction(fraction, octave=2):
    while fraction >= 1:
        fraction = fraction / octave
    while fraction < 1:
        fraction = fraction * octave
    return fraction


def to_octave_reduced_fractions(fractions, octave_size=2):
    reduced_fractions = set()
    for fraction in fractions:
        reduced_fractions.add(to_octave_reduced_fraction(fraction, octave_size))
    return reduced_fractions


def to_octave_reduced_fraction(fraction, octave_size=2):
    while fraction > octave_size:
        fraction /= octave_size
    while fraction < octave_size:
        fraction *= octave_size
    fraction /= octave_size
    return fraction


def generate_matching_chords(chord: set[Fraction, ...]) -> set[tuple[Fraction, ...], ...]:
    # chord 'a' matches chord 'b' if they share a fraction relative to some fundamental
    all_matching_chords = set()
    inverted_fractions = get_inverted_fractions(*chord)
    for inverted_fraction in inverted_fractions:
        for fraction in chord:
            scaling = inverted_fraction * fraction
            all_matching_chords.add(tuple(note * scaling for note in chord))
    return all_matching_chords


def get_overtones_undertones_above_fundamental(number_of_overtones):
    all_tones = set()
    for overtone in range(1, number_of_overtones + 1):
        for undertone in range(1, overtone + 1):
            tone = get_octave_reduction(Fraction(overtone, undertone))
            all_tones.add(tone)
    return all_tones


def get_dissonance_for_fractions(fractions: set[Fraction, ...]) -> Fraction:
    # The dissonance of a set of fractions is the fraction of unique multiples to the fractions' lcm.
    # e.g. fraction 1, 2/3 and 4/5 have lcm 4, with 4 + 6 + 5 - 2 unique multiples
    # The dissonance is thus 'unique multiples'/'lcm' = 12/4 = 3
    fractions_lcm = get_lcm_for_fractions(fractions)
    all_multiples = set()
    for fraction in fractions:
        i = 1
        fraction_multiple = fraction * i
        while fraction_multiple <= fractions_lcm:
            all_multiples.add(fraction_multiple)
            i += 1
            fraction_multiple = fraction * i
    return Fraction(len(all_multiples), fractions_lcm)


def get_consonance_for_fractions(fractions: set[Fraction, ...]) -> Fraction:
    # The consonance of a set of fractions is the number of overlapping multiples per lcm.
    # e.g. fraction 1, 2/3 and 4/5 have lcm 4, with 3 overlaps at 4 and 1 at 2
    # The consonance is thus 1
    fractions_lcm = get_lcm_for_fractions(fractions)
    all_multiples = defaultdict(set)
    for fraction in fractions:
        i = 1
        fraction_multiple = fraction * i
        while fraction_multiple <= fractions_lcm:
            all_multiples[fraction].add(fraction_multiple)
            i += 1
            fraction_multiple = fraction * i

    overlaps = 0
    for comb in combinations(all_multiples.values(), 2):
        overlaps += len(comb[0] & comb[1])

    return Fraction(overlaps, fractions_lcm)


def get_chords_in_scale(scale, min_chord_size=3, max_chord_size=7):
    all_matching_chords = []
    all_non_matching_chords = []

    for chord_size in range(min_chord_size, max_chord_size + 1):
        index_combos = combinations(range(len(scale)), chord_size)
        for index_combo in index_combos:
            for index in index_combo:
                if not scale[index]:
                    all_non_matching_chords.append(index_combo)
                    break
            else:
                all_matching_chords.append(index_combo)
    return all_matching_chords, all_non_matching_chords


def get_chord_matches(chord: tuple, scale: list[int, ...]) -> tuple[str, ...]:
    # print(get_chord_matches((0,4,7),Scale.MAJOR.value))
    # --> ('C', 'F', 'G')
    current_scale = deque(scale)
    matches = []

    for i in range(len(scale)):
        for index in chord:
            if not current_scale[index]:
                break
        else:
            matches.append(NOTE_NAMES[i])
        current_scale.rotate()

    return tuple(matches)


def get_chord_scales(chord: tuple):
    """
    scales = get_chord_scales((0, 4, 7))
    for scale in scales:
        print(scale)
    -->('BLUES', 'C')
    -->('BLUES', 'A')
    -->('HARMONIC_MAJOR', 'G')
    -->...
    """
    matches = []
    for scale in Scale:
        chord_matches = get_chord_matches((0, 3, 7), scale.value)
        for match in chord_matches:
            matches.append((scale.name, match))

    return sorted(matches, key=lambda x: x[0])
