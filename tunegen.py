import os
from collections import defaultdict
from fractions import Fraction
from mautils import get_closest_scientific_pitch

CENT = 2 ** (1 / 1200)
C_4 = 261.6256
TWELVE_TET_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
TWELVE_TET = {TWELVE_TET_NAMES[i]: 2 ** (i / 12) for i in range(12)}
JK_TUNINGS_PATH = r'C:\Users\Johan\AppData\Roaming\REAPER\Data\jk_tunings' + '\\'


# MIDI note 48 is C3
def to_file(frequencies, file_name, first_midi_note=48, path=''):
    frequencies = notes_to_midi(first_midi_note, frequencies)

    file_path = path + file_name
    with open(file_path, 'w') as file:
        for frequency in frequencies[:-1]:
            file.write(str(frequency) + '\n')
        file.write(str(frequencies[-1]))

    if path == '':
        path = os.getcwd()
    path = path.replace('\\', '/')
    file_path = file_path.replace('\\', '/')
    print('Wrote to file "' + file_name + '" file:///' + file_path + ' in folder file:///' + path)


def notes_to_midi(first_midi_note, notes):
    notes = sorted(notes)
    first_non_zero_frequency = next((i for i, x in enumerate(notes) if x), None)
    notes = notes[first_non_zero_frequency:]
    notes = [0] * first_midi_note + notes
    notes += [0] * (128 - len(notes))
    notes = notes[:128]
    return notes


def get_equal_temperament_for_midi(octave_frequency_ratio, keys_per_octave):
    keys_per_octave = float(keys_per_octave)
    midi_note_69 = 440.0
    all_128_midi_notes = [None] * 128
    for i in range(128):
        all_128_midi_notes[i] = midi_note_69 * (octave_frequency_ratio ** ((float(i) - 69.0) / keys_per_octave))
    return all_128_midi_notes


def get_recursive_harmonic_series(fundamental_frequency, number_of_harmonics, recursions):
    harmonic_series = [None] * (number_of_harmonics + (number_of_harmonics - 1) * recursions)
    current_fundamental = fundamental_frequency
    for recursion in range(recursions + 1):
        for key in range(number_of_harmonics):
            harmonic_series[key + (number_of_harmonics - 1) * recursion] = current_fundamental * (key + 1)
        current_fundamental *= number_of_harmonics
    return harmonic_series


def get_harmonic_series(fundamental_frequency, number_of_harmonics):
    harmonic_series = [None] * number_of_harmonics
    current_fundamental = fundamental_frequency
    for key in range(number_of_harmonics):
        harmonic_series[key] = current_fundamental * (key + 1)
    return harmonic_series


def to_octave_reduced_fraction(numerator, octave_size=2):
    octave = 1
    while octave < numerator:
        octave *= octave_size
    if octave != 1:
        octave //= octave_size
    return Fraction(numerator, octave)


def get_octaves(fractions, fundamental, octave_size=2, c_octave=2, octaves=5):
    fractions = sorted(fractions)

    def __fractions_to_keys(keys, keys_in_octave=12):
        midi_notes = [0] * 128
        first_note_index = 12 * (1 + c_octave)
        for octave in range(0, octaves):
            for i, fraction in enumerate(fractions):
                note_index = first_note_index + keys[i]
                if note_index > 127:
                    return midi_notes
                midi_notes[note_index] = fundamental * fraction * (
                        octave_size ** octave)
            first_note_index = first_note_index + keys_in_octave
        return midi_notes

    if len(fractions) < 8:
        white_keys = [0, 2, 4, 5, 7, 9, 11]
        return __fractions_to_keys(white_keys)
    elif len(fractions) < 13:
        all_keys = list(range(12))
        return __fractions_to_keys(all_keys)
    else:
        return __fractions_to_keys(list(range(len(fractions))), len(fractions))


def make_24_et():
    base_frequency = 261.63
    my_keyboard = [base_frequency * 2 ** (i / 24) for i in range(64)]
    to_file(my_keyboard, "24_ET_3_octaves.txt", path=JK_TUNINGS_PATH)


def get_undertones_for_harmonics(number_of_undertones, number_of_harmonics) -> set[Fraction, ...]:
    all_tones = set()

    for harmonic in range(1, number_of_harmonics + 1):
        all_tones = all_tones | get_undertones_for_harmonic(harmonic, number_of_undertones)

    return all_tones


def get_fractions_with_denominator(denominator, max_numerator):
    return {Fraction(i, denominator) for i in range(1, max_numerator) if
            Fraction(i, denominator).denominator == denominator}


def get_fractions(max_numerator, max_denominator):
    fractions = set()
    for i in range(1, max_numerator + 1):
        for j in range(1, max_denominator + 1):
            fractions.add(Fraction(i, j))
    return fractions


def get_undertones_for_harmonic(harmonic, number_of_undertones):
    return {Fraction(harmonic, i) for i in range(1, number_of_undertones + 1)}


def get_notes_by_remainder_dissonance(max_numerator, max_denominator):
    notes = defaultdict(set)
    for numerator in range(1, max_numerator + 1):
        for denominator in range(1, max_denominator + 1):
            fraction = Fraction(numerator, denominator)
            notes[(fraction.numerator - 1) % fraction.denominator].add(fraction)
    return notes


def fractions_dict_to_set(fractions_dict):
    return {fraction for fractions in fractions_dict.values() for fraction in fractions}


def sorted_by_consonance(unique_notes):
    # Note that instead of 1/k, k is used for sorting by consonance to ease indexing
    consonance_sorted = defaultdict(list)
    for fractions in unique_notes.values():
        for fraction in fractions:
            consonance_sorted[fraction.my_denominator].append(fraction)
    return consonance_sorted


def fraction_to_closest_scientific_pitch(fraction, fundamental_frequency=C_4):
    fraction_as_pitch = fraction * fundamental_frequency
    smallest_diff = 100000
    name_number = '-'
    for n in range(-48, 60):
        scientific_pitch = C_4 * 2 ** (n / 12)
        current_diff = abs(fraction_as_pitch - scientific_pitch)
        if current_diff < smallest_diff:
            smallest_diff = current_diff
            name_number = n
    pitch_name = TWELVE_TET_NAMES[name_number % 12] + f'{4 + name_number // 12}'
    return fraction, pitch_name, smallest_diff / (C_4 * 2 ** (name_number / 12)), smallest_diff


def get_fractions_with_absolute_dissonance(dissonance, max_consonance):
    fractions = set()
    for a in range(1, max_consonance + 1):
        fractions.add(Fraction(dissonance + 1, a))
    return fractions


def get_fractions_with_relative_dissonance(relative_dissonance, max_consonance):
    fractions = set()
    for a in range(1, max_consonance + 1):
        fractions.add(Fraction(relative_dissonance * a + 1, a))
    return fractions


def fractions_to_frequencies(fractions, fundamental_tone=261.63):
    return [fundamental_tone * fraction for fraction in fractions]


def make_keyboard_from_fractions(fractions, name="custom_fractions"):
    my_keyboard, first_midi_note = get_centered_keyboard(fractions)
    to_file(fractions_to_frequencies(my_keyboard), f"{name}.txt", first_midi_note=first_midi_note, path=JK_TUNINGS_PATH)

    notes_as_midi = notes_to_midi(first_midi_note, list(my_keyboard))
    print("Keyboard layout:")
    for i in range(10):
        octave_keys = notes_as_midi[12 * i:12 * (i + 1)]
        if any(octave_keys):
            print(f"octave {i - 1}")
            print([get_closest_scientific_pitch(octave_key)[:3] for octave_key in octave_keys if octave_key != 0])


def make_keyboard_with_absolute_dissonance(absolute_dissonance, max_consonance):
    fractions = get_fractions_with_absolute_dissonance(absolute_dissonance, max_consonance)
    fractions.add(Fraction(1, 1))
    my_keyboard, first_midi_note = get_centered_keyboard(fractions)
    to_file(fractions_to_frequencies(my_keyboard), f"abs_diss_{absolute_dissonance}_max_con_{max_consonance}.txt",
            first_midi_note=first_midi_note,
            path=JK_TUNINGS_PATH)
    print(sorted(fractions))


def make_keyboard_with_relative_dissonance(relative_dissonance, max_consonance):
    fractions = get_fractions_with_relative_dissonance(relative_dissonance, max_consonance)
    fractions.add(Fraction(1, 1))
    my_keyboard, first_midi_note = get_centered_keyboard(fractions)
    to_file(fractions_to_frequencies(my_keyboard), f"rel_diss_{relative_dissonance}_max_con_{max_consonance}.txt",
            first_midi_note=first_midi_note, path=JK_TUNINGS_PATH)
    print(sorted(fractions))


def get_fractions_with_remainder_dissonance(remainder_dissonance, max_consonance):
    fractions = set()
    for a in range(1, max_consonance + 1):
        for b in range(1, 3 * a + 1):
            if (b - 1) % a == remainder_dissonance:
                fractions.add(Fraction(b, a))
    return fractions


def make_keyboard_with_remainder_dissonance(remainder_dissonance, max_consonance):
    fractions = get_fractions_with_remainder_dissonance(remainder_dissonance, max_consonance)
    fractions.add(Fraction(1, 1))
    my_keyboard, first_midi_note = get_centered_keyboard(fractions)
    to_file(fractions_to_frequencies(my_keyboard), f"rem_diss_{remainder_dissonance}_max_con_{max_consonance}.txt",
            first_midi_note=first_midi_note, path=JK_TUNINGS_PATH)
    print(sorted(fractions))


def make_keyboard_from_fractions_with_denominator(denominator, number_of_fractions):
    fractions = get_fractions_with_denominator(denominator, number_of_fractions)
    fractions.add(Fraction(1, 1))
    my_keyboard, first_midi_note = get_centered_keyboard(fractions)
    to_file(fractions_to_frequencies(my_keyboard), f"denominator_{denominator}_notes_{number_of_fractions}.txt",
            first_midi_note=first_midi_note, path=JK_TUNINGS_PATH)
    print(sorted(fractions))


def make_keyboard_from_harmonic_undertones(max_harmonic, max_undertones):
    fractions = set()
    for harmonic in range(1, max_harmonic + 1):
        fractions.update(get_undertones_for_harmonic(harmonic, max_undertones))
    fractions.add(Fraction(1, 1))
    my_keyboard, first_midi_note = get_centered_keyboard(fractions)
    to_file(fractions_to_frequencies(my_keyboard), f"max_harm_{max_harmonic}_undertones_{max_undertones}.txt",
            first_midi_note=first_midi_note, path=JK_TUNINGS_PATH)
    print(sorted(fractions))


def get_centered_keyboard(fractions: set, c_octave=4, keyboard_size=60):
    fractions.add(Fraction(1, 1))
    fractions = sorted(list(fractions))
    index_of_fundamental = fractions.index(Fraction(1, 1))

    if len(fractions) > keyboard_size:
        fractions = fractions[index_of_fundamental - keyboard_size // 2:index_of_fundamental + keyboard_size // 2]
        index_of_fundamental = fractions.index(Fraction(1, 1))

    lower_frequencies = [fraction for fraction in fractions[:index_of_fundamental]]
    higher_frequencies = [fraction for fraction in fractions[index_of_fundamental:]]
    center_midi = 12 * (c_octave + 1)
    return lower_frequencies + higher_frequencies, center_midi - len(lower_frequencies)


def sort_by_scientific_notation(data_list):
    def __sort_data():
        for name in TWELVE_TET_NAMES:
            for number in range(12):
                if name + str(number) in data[1]:
                    sorted_data[name].append(data)
                    return

    sorted_data = defaultdict(list)
    for data in data_list:
        __sort_data()

    return {name: sorted(sorted_data[name], key=lambda x: x[2]) for name in TWELVE_TET_NAMES}
