import math
from collections import defaultdict, deque
from fractions import Fraction
from itertools import combinations, product
from math import gcd, lcm
from sympy.ntheory import factorint
from enum import Enum


class Scale(Enum):
    ARABIAN = [1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0]
    BLUES = [1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0]  # blues minor
    BYZANTINE = [1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1]
    # BLUES_EXTENDED is more related to A# than C?
    BLUES_EXTENDED = [1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0]
    HARMONIC_MAJOR = [1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1]  # C + Fm + G
    HARMONIC_MINOR = [1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1]  # Cm + Fm + G
    HARMONIC_MINOR_EXTENDED = [1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1]
    MAJOR = [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1]  # C + F + G
    MELODIC_MINOR = [1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1]  # Cm + F + G
    PENTATONIC = [1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0]  # C Major
    HEXATONIC = [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0]  # C Major
    OCTACTONIC = [1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0]
    # C Natural scale, based on dissonance curve minima for harmonic overtones.
    # See the book "Tuning, Timbre, Spectrum, Scale".
    NATURAL_SCALE = [1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0]


MAJOR_SCALE_FRACTIONS = {Fraction(1), Fraction(9, 8), Fraction(5, 4), Fraction(4, 3), Fraction(3, 2), Fraction(5, 3),
                         Fraction(15, 8)}
NOTE_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']


def identify_scale(scale_to_identify: list) -> (list, str):
    # print(identify_scale([1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1]))
    # --> ('MAJOR', 'C')
    scale_to_identify = deque(scale_to_identify)
    fundamental = 0
    for i in range(len(scale_to_identify)):
        for scale in Scale:
            if list(scale_to_identify) == scale.value:
                return scale.name, NOTE_NAMES[fundamental]
        scale_to_identify.rotate()
        fundamental -= 1


def parse_scales(path):
    """
    for my_scale, my_comment in parse_scales(r'D:\Projects\Code\MusicAnalysis\scale.txt'):
        if my_scale is None:
            print(my_comment)
        else:
            print(identify_scale(my_scale), my_comment)
    """
    with open(path, 'r') as file:
        output = []

        for line in file:
            line = line.strip('\n').split(' ')
            comment = ''

            # ignore comments after //
            for index, note in enumerate(line):
                if '//' in note:
                    comment = " ".join(line[index + 1:])
                    line = line[:index]
                    break

            if len(line) == 0:
                output.append((None, comment))
                continue

            # format 'b' to '#'
            for index, note in enumerate(line):
                if 'b' in note:
                    line[index] = NOTE_NAMES[NOTE_NAMES.index(note[0]) - 1]

            # create scale
            scale = [0] * 12
            for note in line:
                if note in NOTE_NAMES:
                    scale[NOTE_NAMES.index(note)] = 1

            output.append((scale, comment))
    return output


def get_sensory_dissonance(frequency_series, beating_bandwidth=20.0):
    """
    takes a list of frequencies and calculates their sensory dissonance
    :return: total number of dissonant intervals
    """
    dissonances = []
    for frequency_serie_index in range(len(frequency_series)):
        frequency_serie = frequency_series[frequency_serie_index]
        for frequency_index, frequency in enumerate(frequency_serie):
            for other_frequency_serie_index in range(frequency_serie_index + 1, len(frequency_series)):
                other_frequency_serie = frequency_series[other_frequency_serie_index]
                for other_frequency_index, other_frequency in enumerate(other_frequency_serie):
                    interval_abs = abs(frequency - other_frequency)
                    if beating_bandwidth < interval_abs < get_ERB(frequency) - (beating_bandwidth / 2.0):
                        dissonances.append((frequency, other_frequency, frequency_index, other_frequency_index,
                                            frequency_serie_index, other_frequency_serie_index))
    return dissonances


def get_ERB(hz):
    # https://en.wikipedia.org/wiki/Equivalent_rectangular_bandwidth
    # returns ERB in hz
    khz = hz / 1000.0
    return 6.23 * khz * khz + 93.39 * khz + 28.53


def get_harmonic_series(fundamentals, series_length):
    all_series = []
    for fundamental in fundamentals:
        all_series.append([i * fundamental for i in range(1, series_length + 1)])
    return all_series


def equal_under_rotation(scale_1, scale_2) -> bool:
    scale_1 = deque(scale_1)
    scale_2 = deque(scale_2)
    for i in range(len(scale_1)):
        if scale_1 == scale_2:
            return True
        scale_1.rotate()
    return False


def get_all_12_tet_scales(scale_size: int, max_consecutive_notes: int) -> list[list[int, ...]]:
    """
    for my_scale in get_all_12_tet_scales(7, 3):
        print(my_scale)
    --> [1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1]
    --> ...
    """
    all_scale_configurations = []
    for prod in product([0, 1], repeat=11):
        prod = deque([1]) + deque(prod)
        all_scale_configurations.append(prod)

    allowed_scales = []
    for current_scale in all_scale_configurations:
        if sum(current_scale) != scale_size:
            continue

        if get_max_consecutive_notes(current_scale) > max_consecutive_notes:
            continue
        else:
            allowed_scales.append(current_scale)

    unique_scales = []
    for current_scale in allowed_scales:
        for unique_scale in unique_scales:
            if equal_under_rotation(current_scale, unique_scale):
                break
        else:
            unique_scales.append(list(current_scale))

    return unique_scales


def get_number_of_semitone_intervals(notes):
    """
    notes = [1,0,1,1]
    print(get_number_of_semitone_intervals(notes))
    --> 2
    """
    semitone_intervals = 0
    for i in range(len(notes)):
        if notes[i] == notes[(i + 1) % len(notes)] == 1:
            semitone_intervals += 1
    return semitone_intervals


def get_max_consecutive_notes(notes):
    """
    notes = [1,1,1,0,1]
    print(get_max_consecutive_notes(notes))
    --> 4
    """
    max_consecutive_notes = 0
    current_consecutive_notes = 0
    for i in range(2 * len(notes)):
        if notes[i % len(notes)] == 1:
            current_consecutive_notes += 1
            if current_consecutive_notes > max_consecutive_notes:
                max_consecutive_notes = current_consecutive_notes
        else:
            current_consecutive_notes = 0
    return max_consecutive_notes


def get_closest_scientific_pitch(fraction):
    a_4 = 440.0
    twelve_tet_names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    fraction_as_pitch = fraction * a_4
    smallest_diff = 100000
    name_number = '-'
    for n in range(-39, 69):
        scientific_pitch = a_4 * 2 ** (n / 12)
        current_diff = abs(fraction_as_pitch - scientific_pitch)
        if current_diff < smallest_diff:
            smallest_diff = current_diff
            name_number = n
    pitch_name = twelve_tet_names[name_number % 12] + f'{4 + name_number // 12}'
    closest_pitch = a_4 * 2 ** (name_number / 12)
    return fraction, pitch_name, 1200 * math.log2(
        fraction_as_pitch / closest_pitch), smallest_diff / closest_pitch, smallest_diff


def get_midi_pitches():
    concert_pitch = 440
    interval_ratio = 2 ** (1 / 12.0)
    midi_pitches = []
    current_pitch = concert_pitch
    for i in range(1, 70):
        current_pitch /= interval_ratio
        midi_pitches.append(current_pitch)
    current_pitch = concert_pitch
    midi_pitches.sort()
    for i in range(59):
        midi_pitches.append(current_pitch)
        current_pitch *= interval_ratio
    return midi_pitches


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
    """
    major_chord_C = [Fraction(1), Fraction(5, 4), Fraction(3, 2)]
    for value in get_possible_lcm_configurations_for_fractions(*major_chord_C):
        print(value)
    -->
    (Fraction(1, 1), Fraction(5, 4), Fraction(3, 2))
    (Fraction(15, 1), Fraction(12, 1), Fraction(10, 1))
    ({Fraction(1, 1), Fraction(5, 4), Fraction(3, 2)}, {Fraction(1, 1), Fraction(6, 5), Fraction(4, 5)}, {Fraction(1, 1), Fraction(2, 3), Fraction(5, 6)})
    """
    fractions = sorted(set(fractions))
    lcm_configurations = []
    lcm_values = []

    for reference_fraction in fractions:
        lcm_configurations.append({fraction / reference_fraction for fraction in fractions})
        lcm_values.append(get_lcm_for_fractions(lcm_configurations[-1]))

    return tuple(fractions), tuple(lcm_values), tuple(lcm_configurations)


def get_lcd_for_fractions(fractions: set[Fraction, ...], ignore_zero=True) -> Fraction:
    if ignore_zero:
        fractions = [fraction for fraction in fractions if fraction != 0]
    fractions_gcd = gcd(*[fraction.numerator for fraction in fractions])
    fractions_lcd = Fraction(lcm(*[fraction.denominator for fraction in fractions]), fractions_gcd)
    return fractions_lcd


def get_possible_lcd_configurations_for_fractions(*fractions: Fraction):
    """
    major_chord_C = [Fraction(1), Fraction(5, 4), Fraction(3, 2)]
    for value in get_possible_lcd_configurations_for_fractions(*major_chord_C):
        print(value)
    -->
    (Fraction(1, 1), Fraction(5, 4), Fraction(3, 2)) <-- original fractions
    (Fraction(15, 1), Fraction(12, 1), Fraction(10, 1)) <-- LCD values for...
    ({Fraction(1, 1), Fraction(5, 4), Fraction(3, 2)}, {Fraction(1, 1), Fraction(6, 5), Fraction(4, 5)}, {Fraction(1, 1), Fraction(2, 3), Fraction(5, 6)}) <-- ...the fraction permutations
    """
    fractions = sorted(set(fractions))
    lcd_configurations = []
    lcd_values = []

    for reference_fraction in fractions:
        lcd_configurations.append({fraction / reference_fraction for fraction in fractions})
        lcd_values.append(get_lcd_for_fractions(lcd_configurations[-1]))

    return tuple(fractions), tuple(lcd_values), tuple(lcd_configurations)


def get_inverted_fractions(*fractions: Fraction) -> list[Fraction, ...]:
    return [Fraction(fraction.denominator, fraction.numerator) for fraction in fractions]


# TODO return with chord[0] = 1, all scales distinct under rotation
def get_all_12_tet_chords():
    # Chord defined as any set of octave reduced notes in root position which does not contain a semitone intervall
    all_key_combinations = product((0, 1), repeat=9)
    relevant_key_combinations = []
    for comb in all_key_combinations:
        for i in range(len(comb) - 1):
            if comb[i] == comb[i + 1] == 1:
                break
        else:
            relevant_key_combinations.append((1, 0, *comb, 0))

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
        chord_matches = get_chord_matches(chord, scale.value)
        for match in chord_matches:
            matches.append((scale.name, match))

    return sorted(matches, key=lambda x: x[0])
