import math
from collections import defaultdict
from fractions import Fraction
import matplotlib.pyplot as plt
from math import lcm
from math import gcd
from itertools import combinations
from matplotlib.patches import Rectangle

C_4 = 261.6256
TWELVE_TET_NAMES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']


def get_closest_scientific_pitch(fraction):
    fraction_as_pitch = fraction * C_4
    smallest_diff = 100000
    name_number = '-'
    for n in range(-48, 60):
        scientific_pitch = C_4 * 2 ** (n / 12)
        current_diff = abs(fraction_as_pitch - scientific_pitch)
        if current_diff < smallest_diff:
            smallest_diff = current_diff
            name_number = n
    pitch_name = TWELVE_TET_NAMES[name_number % 12] + f'{4 + name_number // 12}'
    closest_pitch = C_4 * 2 ** (name_number / 12)
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


def get_all_notes_with_dissonance_less_than(max_dissonance):
    all_notes = defaultdict(list)
    for dissonance in range(2, max_dissonance + 1):
        for j in range(1, dissonance):
            k = dissonance - j
            fraction = Fraction(j, k)
            if fraction not in all_notes[dissonance]:
                all_notes[dissonance].append(fraction)
    for notes in all_notes.values():
        notes.sort()
    return all_notes


def get_unique_notes_with_dissonance_less_than(max_dissonance):
    unique_notes = {fraction for fractions in get_all_notes_with_dissonance_less_than(max_dissonance).values() for
                    fraction in fractions}
    return sorted(list(unique_notes))


def get_colors_by_consonance(unique_fractions, max_dissonance):
    rgb = [(0, 0, 0), (0, 0, 1), (0, 0.5, 0.5), (0, 1, 0), (0.5, 0.5, 0), (1, 0, 0)] + [(0, 0, 0)] * max_dissonance
    return [rgb[fraction.denominator] for fraction in unique_fractions]


def get_heights_by_dissonance(unique_fractions):
    return [fraction.numerator + fraction.denominator - 1 for fraction in unique_fractions]


def get_widths_by_fractions(unique_fractions):
    widths = [(unique_fractions[i + 1] - unique_fractions[i]) * 0.8 / 2 for i in range(len(unique_fractions) - 1)]
    widths = widths + [widths[-1]]
    return widths


def get_undertones_for_harmonic(harmonic, number_of_undertones):
    return [Fraction(harmonic, i) for i in range(1, number_of_undertones + 1)]


def get_overtones_for_fraction(fraction: Fraction, overtones: int = 0) -> list[Fraction, ...]:
    return [fraction + i * fraction for i in range(overtones + 1)]


def get_overtones_for_fractions(overtones: int = 0, *fractions: Fraction) -> set[Fraction, ...]:
    return {overtone for fraction in fractions for overtone in get_overtones_for_fraction(fraction, overtones)}


def get_octave_reduction(fraction, octave=2):
    while fraction >= 1:
        fraction = fraction / octave
    while fraction < 1:
        fraction = fraction * octave
    return fraction


def get_octave_pair(fraction, octave=2):
    fraction = get_octave_reduction(fraction, octave)
    octave_pair = Fraction(fraction.denominator, fraction.numerator)
    return fraction, get_octave_reduction(octave_pair, octave)


def get_extended_octave_pair(fraction, octaves):
    extended_octave_pair = set()
    for octave in octaves:
        extended_octave_pair.update(get_octave_pair(fraction, octave))
    return extended_octave_pair


def plot_notes_with_dissonance_less_than(max_dissonance):
    # noinspection PyTypeChecker
    fig, axs = plt.subplots(1, 2, sharey=True)
    ratios = get_unique_notes_with_dissonance_less_than(max_dissonance)
    ratios_one_and_below = ratios[:ratios.index(Fraction(1, 1)) + 1]
    ratios_one_and_above = ratios[ratios.index(Fraction(1, 1)):]

    ax1 = axs[0]
    ax1.set_xscale('log')
    bar_colors = get_colors_by_consonance(ratios_one_and_below, max_dissonance)
    bar_heights = get_heights_by_dissonance(ratios_one_and_below)
    bar_widths = get_widths_by_fractions(ratios_one_and_below)
    ax1.bar(ratios_one_and_below, bar_heights, width=bar_widths, color=bar_colors,
            **{'edgecolor': 'black', 'linewidth': 0.5})

    ax2 = axs[1]
    bar_colors = get_colors_by_consonance(ratios_one_and_above, max_dissonance)
    bar_heights = get_heights_by_dissonance(ratios_one_and_above)
    bar_widths = get_widths_by_fractions(ratios_one_and_above)
    ax2.bar(ratios_one_and_above, bar_heights, width=bar_widths, color=bar_colors,
            **{'edgecolor': 'black', 'linewidth': 0.5})

    fig.subplots_adjust(wspace=0)

    ax1.set_xlim(float(ratios_one_and_below[0]), float(ratios_one_and_below[-1]))
    ax2.set_xlim(float(ratios_one_and_above[0]), float(ratios_one_and_above[-1]))

    ax1.spines.right.set_visible(False)
    ax2.spines.left.set_visible(False)
    ax2.tick_params(which='both', left=False)
    plt.setp(ax2.get_xticklabels()[0], visible=False)

    plt.show()


def plot_wavelength_multiples_for_fraction_sets(fraction_sets: list[set[Fraction, ...], ...], lcm_plot=True,
                                                title=None, label_style='None', plot_style='None'):
    fig, axs = plt.subplots(len(fraction_sets), 1, layout='constrained', squeeze=False)
    for i in range(len(fraction_sets)):
        plot_wavelength_multiples_for_fractions(axs[i][0], fraction_sets[i], lcm_plot, title, label_style, plot_style)


def plot_wavelength_multiples_for_fractions(ax, fractions, lcm_plot, title, label_style, plot_style):
    wavelengths = sorted([Fraction(fraction.denominator, fraction.numerator) for fraction in fractions], reverse=True)
    global_wavelength_lcm = get_lcm_for_fractions(*wavelengths)
    all_multiples = set()
    multiples_for_wavelength = {}

    for wavelength in wavelengths:
        if lcm_plot:
            multiples = int(global_wavelength_lcm * wavelength.denominator / wavelength.numerator)
        else:
            multiples = wavelength.denominator

        current_wavelengths = [wavelength * i for i in range(1, multiples + 1)]
        multiples_for_wavelength[wavelength] = current_wavelengths
        all_multiples.update(current_wavelengths)

    if plot_style == 'compact':
        new_multiples_for_wavelength = {}
        old_multiples_for_wavelength = {key: value for key, value in multiples_for_wavelength.items()}
        for wavelength, multiples in multiples_for_wavelength.items():
            del old_multiples_for_wavelength[wavelength]
            for old_multiples in old_multiples_for_wavelength.values():
                if set(multiples).issubset(set(old_multiples)):
                    break
            else:
                new_multiples_for_wavelength[wavelength] = multiples

        multiples_for_wavelength = new_multiples_for_wavelength

    for wavelength in multiples_for_wavelength.keys():
        current_wavelengths = multiples_for_wavelength[wavelength]
        bar_heights = [wavelength] * len(current_wavelengths)
        bar_widths = [-wavelength] * len(current_wavelengths)
        ax.bar(current_wavelengths, bar_heights, bar_widths, align='edge', **{'edgecolor': 'black', 'linewidth': 0.5})

    if label_style == 'lcm':
        lcm_for_combinations = get_lcm_for_combinations(set(wavelengths))
        lcm_multiples = set()
        for lcm_for_combination in lcm_for_combinations:
            multiples = int(global_wavelength_lcm * lcm_for_combination.denominator / lcm_for_combination.numerator)
            lcm_multiples.update([lcm_for_combination * i for i in range(1, multiples + 1)])

        lcm_multiples.add(Fraction(0))
        x_ticks = list(lcm_multiples)
    else:
        x_ticks = sorted(all_multiples) + [0]
    y_ticks = sorted(wavelengths) + [Fraction(0), Fraction(1)]

    pitch_classes = set()
    for tick in sorted(x_ticks, key=lambda fraction: (fraction.denominator, - fraction.numerator)):
        try:
            scientific_pitch = get_closest_scientific_pitch(Fraction(tick.denominator, tick.numerator))
            pitch_classes.add(scientific_pitch[1][:-1])
            print(scientific_pitch)
        except ZeroDivisionError:
            pass
    print(len(pitch_classes), sorted(pitch_classes))
    print()
    for wavelength in sorted(wavelengths):
        try:
            scientific_pitch = get_closest_scientific_pitch(Fraction(wavelength.denominator, wavelength.numerator))
            pitch_classes.add(scientific_pitch[1][:-1])
            print(scientific_pitch)
        except ZeroDivisionError:
            pass
    print(len(pitch_classes), sorted(pitch_classes))
    print('---')
    print()

    ax.set_xticks(x_ticks, x_ticks)
    ax.set_yticks(y_ticks, y_ticks)

    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Original Wavelength')
    ax.margins(0)
    ax.set_title(title)


def temp():
    A_sharp_6_major = {Fraction(7), Fraction(7 * 5, 4), Fraction(7 * 3, 2)}
    C_4_major = {Fraction(1), Fraction(5, 4), Fraction(3, 2)}
    C_5_major = {Fraction(2), Fraction(2 * 5, 4), Fraction(2 * 3, 2)}
    C_6_major = {Fraction(4), Fraction(4 * 5, 4), Fraction(4 * 3, 2)}
    C_4_minor = {Fraction(1), Fraction(6, 5), Fraction(3, 2)}
    C_sharp_4_major = {Fraction(16, 15), Fraction(16 * 5, 15 * 4), Fraction(16 * 3, 15 * 2)}
    D_4_major = {Fraction(9, 8), Fraction(9 * 5, 8 * 4), Fraction(9 * 3, 8 * 2)}
    E_6_major = {Fraction(5), Fraction(5 * 5, 4), Fraction(5 * 3, 2)}
    F_4_major = {Fraction(4, 3), Fraction(4 * 5, 3 * 4), Fraction(4 * 3, 3 * 2)}
    F_3_major = {Fraction(2, 3), Fraction(2 * 5, 3 * 4), Fraction(2 * 3, 3 * 2)}
    G_4_major = {Fraction(3, 2), Fraction(15, 8), Fraction(9, 4)}
    G_5_major = {Fraction(3), Fraction(3 * 5, 4), Fraction(3 * 3, 2)}
    G_6_major = {Fraction(6), Fraction(6 * 5, 4), Fraction(6 * 3, 2)}
    Bb_4_major = {Fraction(7, 4), Fraction(7 * 5, 4 * 4), Fraction(7 * 3, 4 * 2)}
    Bb_6_major = {Fraction(7, 1), Fraction(7 * 5, 1 * 4), Fraction(7 * 3, 1 * 2)}

    my_fractions = [C_4_major, C_4_major | C_5_major, C_4_major | C_6_major]
    """
    my_fractions = [
        # C_4_major,
        # C_4_major | C_5_major,
        C_4_major | C_5_major | G_5_major,
        # C_4_major | C_5_major | G_5_major | C_6_major,
        # C_4_major | C_5_major | G_5_major | C_6_major | E_6_major,
        # C_4_major | C_5_major | G_5_major | C_6_major | E_6_major | G_6_major,
        # C_4_major | C_5_major | G_5_major | C_6_major | E_6_major | G_6_major | A_sharp_6_major
    ]

    my_inverted_fractions = [{Fraction(fraction.denominator, fraction.numerator) for fraction in fractions} for
                             fractions in my_fractions]
    """
    # plot_wavelength_multiples_for_fraction_sets(my_fractions)
    max_harmonic = 14
    my_fractions = [
        # get_overtones_for_fractions(19, Fraction(1)),
        # get_overtones_for_fractions(8, Fraction(1)),
        get_overtones_for_fractions(8, Fraction(1)) | get_overtones_for_fractions(11, Fraction(2, 3)),
        get_overtones_for_fractions(11, Fraction(1)) | get_overtones_for_fractions(8, Fraction(3, 2)),

        # get_overtones_for_fractions(5, Fraction(1)) | get_overtones_for_fractions(8, Fraction(1, 2)),
        # get_overtones_for_fractions(8, Fraction(1)) | get_overtones_for_fractions(3, Fraction(3, 2)),
        # get_overtones_for_fractions(11, Fraction(1)) | get_overtones_for_fractions(9, Fraction(6, 5)),
        # get_overtones_for_fractions(8, Fraction(1)) | get_overtones_for_fractions(8, Fraction(2, 3)),
        # get_overtones_for_fractions(5, Fraction(1)) | get_overtones_for_fractions(6, Fraction(3, 4)),

        # get_overtones_for_fractions(max_harmonic, Fraction(1), Fraction(3, 2)),
        # get_overtones_for_fractions(max_harmonic, Fraction(1), Fraction(3)),
        # get_overtones_for_fractions(max_harmonic, Fraction(1), Fraction(3, 2), Fraction(3)),
        # get_overtones_for_fractions(15, Fraction(1)) | get_overtones_for_fractions(6, Fraction(3))
    ]
    for fractions in my_fractions:
        plot_wavelength_multiples_for_fraction_sets([fractions], label_style='lcm')
        plot_wavelength_multiples_for_fraction_sets([fractions], plot_style='compact')
        # plot_wavelength_multiples_for_fraction_sets(my_inverted_fractions, lcm_plot=True)


def plot_rectangles(ax, *rectangles: Rectangle):
    for rectangle in rectangles:
        ax.add_patch(rectangle)


def plot_undertones(fractions: list[Fraction], number_of_harmonics: list[int], max_harmonic=9):
    if len(fractions) != len(number_of_harmonics):
        print(f"Number of fractions does not match number of harmonics: {len(fractions)} != {len(number_of_harmonics)}")
        return

    fig, ax = plt.subplots(layout='constrained')
    original_wavelengths = sorted([Fraction(fraction.denominator, fraction.numerator) for fraction in fractions],
                                  reverse=True)
    global_wavelength_lcm = get_lcm_for_fractions(*original_wavelengths)
    number_of_overtones = [harmonic - 1 for harmonic in number_of_harmonics]

    # init wavelengths_by_denominator
    wavelengths_by_denominator = defaultdict(list)
    for i in range(len(fractions)):
        for overtone in get_overtones_for_fraction(fractions[i], number_of_overtones[i]):
            # overtone numerator is corresponding wavelength's denominator
            if not overtone.numerator > max_harmonic:
                wavelengths_by_denominator[overtone.numerator].append(
                    Fraction(overtone.denominator, overtone.numerator))

    global_alpha = 1 / max([len(wavelengths) for wavelengths in wavelengths_by_denominator.values()])

    # init empty rectangles
    rectangles_by_wavelength_denominator = defaultdict(list)
    current_y = 0
    for denominator in sorted(wavelengths_by_denominator.keys(), reverse=True):
        current_x = 0
        rectangle_size = 1 / denominator
        rectangle_params = {'edgecolor': 'black', 'linewidth': 0.5, 'facecolor': f"C{denominator - 1}"}
        for i in range(int(global_wavelength_lcm * denominator)):
            current_rectangle = Rectangle((current_x, current_y), rectangle_size, rectangle_size - current_y,
                                          **rectangle_params)
            current_rectangle.set_alpha(0)
            rectangles_by_wavelength_denominator[denominator].append(current_rectangle)
            current_x += rectangle_size
        current_y = rectangle_size

    # update wavelength intensity
    for denominator, wavelengths in wavelengths_by_denominator.items():
        for wavelength in wavelengths:
            for multiple in range(1, int(global_wavelength_lcm * wavelength.denominator / wavelength.numerator) + 1):
                current_rectangle = rectangles_by_wavelength_denominator[denominator][
                    wavelength.numerator * multiple - 1]
                current_rectangle.set_alpha(current_rectangle.get_alpha() + global_alpha)

    # draw rectangles
    for rectangles in rectangles_by_wavelength_denominator.values():
        plot_rectangles(ax, *rectangles)

    # set labels etc.
    unique_wavelengths = {wavelength for wavelengths in wavelengths_by_denominator.values() for wavelength in
                          wavelengths}
    lcm_for_combinations = get_lcm_for_combinations(set(unique_wavelengths))
    lcm_multiples = set()
    for lcm_for_combination in lcm_for_combinations:
        multiples = int(global_wavelength_lcm * lcm_for_combination.denominator / lcm_for_combination.numerator)
        lcm_multiples.update([lcm_for_combination * i for i in range(1, multiples + 1)])

    lcm_multiples.add(Fraction(0))
    x_ticks = list(lcm_multiples)

    unique_denominators = {Fraction(1, wavelength.denominator) for wavelength in unique_wavelengths}
    y_ticks = sorted(unique_denominators) + [Fraction(0), Fraction(1)]

    ax.set_xticks(x_ticks, x_ticks)
    ax.set_yticks(y_ticks, y_ticks)

    ax.set_xlabel('In-ear wavelength')
    ax.set_ylabel('Wavelength of harmonic')
    ax.margins(0)
    title = ""
    for fraction in [f"{fraction.numerator}/{fraction.denominator}" for fraction in fractions]:
        title += fraction + ', '
    ax.set_title(title[:-2])


def print_wavelengths_for_harmonics(fraction: Fraction, harmonics: int, max_wavelength_denominator=None):
    if not max_wavelength_denominator:
        max_wavelength_denominator = fraction.numerator * harmonics
    for i in range(1, harmonics + 1):
        current_fraction = fraction * i
        if current_fraction.numerator <= max_wavelength_denominator:
            print(Fraction(current_fraction.denominator, current_fraction.numerator))
    print()


if __name__ == '__main__':

    my_harmonics_max = 20
    my_view_max = 18

    my_chord = [Fraction(1), Fraction(5, 2)]
    plot_undertones(my_chord, [my_harmonics_max] * len(my_chord), my_view_max)

    my_chord = [Fraction(1), Fraction(5, 4)]
    plot_undertones(my_chord, [my_harmonics_max] * len(my_chord), my_view_max)

    my_chord = [Fraction(1), Fraction(5, 3)]
    plot_undertones(my_chord, [my_harmonics_max] * len(my_chord), my_view_max)

    my_chord = [Fraction(1), Fraction(9, 5)]
    plot_undertones(my_chord, [my_harmonics_max] * len(my_chord), my_view_max)

    my_chord = [Fraction(1), Fraction(9, 7)]
    plot_undertones(my_chord, [my_harmonics_max] * len(my_chord), my_view_max)

    my_chord = [Fraction(1), Fraction(9, 8)]
    plot_undertones(my_chord, [my_harmonics_max] * len(my_chord), my_view_max)

    for my_fraction in my_chord[1:]:
        print_wavelengths_for_harmonics(my_fraction, my_view_max, my_view_max)

    plt.show()
