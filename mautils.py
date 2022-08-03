import math


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
