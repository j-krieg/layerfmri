
def unify_color(color):
    # Converts any color oto hex. All different kinds of input are possible:
    #    (52, 152, 219) # Output: '#3498db'
    #    "0.5"          # Output: '#7f7f7f


    if isinstance(color, tuple):
        hex_color = '#{:02x}{:02x}{:02x}'.format(color[0], color[1], color[2])
    elif isinstance(color, str):
        if ',' in color:
            color_values = color.split(',')
            hex_color = '#{:02x}{:02x}{:02x}'.format(int(color_values[0]), int(color_values[1]), int(color_values[2]))
        else:
            hex_color = color
    elif isinstance(color, float):
        color_int = int(color * 255)
        hex_color = '#{:02x}{:02x}{:02x}'.format(color_int, color_int, color_int)
    else:
        return None

    return hex_color


def adjust_color(color, rel_saturation, rel_lightness):
    # adjust saturation and lightness of a color (any input possible, will be unified to hex format)

    # Convert HEX to RGB
    hex_color = unify_color(color) # e.g., '#3498db'
    rgb = tuple(int(hex_color[i:i + 2], 16) for i in (1, 3, 5))

    # RGB to HSL
    h, l, s = colorsys.rgb_to_hls(*[x / 255 for x in rgb])

    # Adjust saturation and lightness
    s = max(0, min(1, s + rel_saturation))
    l = max(0, min(1, l + rel_lightness))

    # HSL to RGB
    r, g, b = [round(x * 255) for x in colorsys.hls_to_rgb(h, l, s)]
    return f'#{r:02x}{g:02x}{b:02x}'




custom_colors = {
    "Curious Blue": (52, 152, 219),
    "Buttercup": (243, 156, 18),
    "Cinnabar": (231, 76, 60),
    "Mountain Meadow": (26, 188, 156),
    "Ripe Lemon": (241, 196, 15),
    "Mariner": (41, 128, 185),
    "Wisteria": (155, 89, 182),
    "Shamrock": (46, 204, 113),
    "Pickled Bluewood": (52, 73, 94),
    "Zest": (230, 126, 34),
    "Tall Poppy": (192, 57, 43),
    "Jungle Green": (39, 174, 96),
    "Burnt Orange": (211, 84, 0),
    "Java": (23, 190, 207),
    "Sunglo": (227, 119, 194),
    "Amaranth": (229, 43, 80),
    "Bittersweet": (254, 111, 94),
    "Brink Pink": (251, 96, 127),
    "Flirt": (162, 0, 109),
    "Jazzberry Jam": (165, 11, 94),
    "Lavender Purple": (150, 123, 182),
    "Lipstick": (171, 5, 99),
    "Monarch": (139, 7, 35),
    "Moody Blue": (127, 118, 211),
    "Mosque": (3, 106, 110),
    "Orient": (1, 94, 133),
    "Paris M": (38, 5, 106),
    "Pelorous": (62, 171, 191),
    "Persian Plum": (112, 28, 28),
    "Radical Red": (255, 53, 94),
    "Sapphire": (47, 81, 158),
    "Seagull": (128, 204, 234),
    "Sunflower": (228, 212, 34),
    "Teal": (0, 128, 128),
    "Terracotta": (226, 114, 91),
    "Wasabi": (120, 138, 37),
    "Wine Berry": (89, 29, 53),
    "West Side": (255, 145, 15),
    "Green Yellow": (173, 255, 47),
    "Mustard": (255, 219, 88),
    "Hacienda": (152, 129, 27),
    "Honeysuckle": (237, 252, 132)
}


# different groups of colours that are most distinguishable between and winthin groups
blues 	= ["Curious Blue", "Java", "Sapphire", "Seagull"]
oranges	= ["West Side", "Mustard", "Honeysuckle", "Hacienda"]
reds	= ["Radical Red", "Cinnabar", "Bittersweet", "Persian Plum"]
greens	= ["Shamrock", "Teal", "Wasabi", "Green Yellow"]
purples	= ["Wisteria", "Sunglo", "Moody Blue", "Jazzberry Jam"]

# selection of 20 well-distinguishable colors for screen and print
col_seq             = [] # list of tupels       e.g., (237, 252, 132)
col_hex             = [] # list of hex strings  e.g., #FF0042
col_names           = [] # list of names        e.g., "Wisteria"
col_dict            = {} # dict wherein names are keys and tuples are values
col_sec_entries     = 20 # e.g., use like that: color = col_seq[n%col_sec_entries]
col_blues           = []
col_oranges         = []
col_reds            = []
col_greens          = []
col_purples         = []

for diff in [0, 1, 2, 3]:
    for n, group in enumerate([blues, oranges, reds, greens, purples]):
        name = group[diff]
        col_tuple = custom_colors[name]

        col_dict[name] = col_tuple
        col_seq.append(col_tuple)
        col_names.append(name)
        col_hex.append(unify_color(col_tuple))

        (([col_blues, col_oranges, col_reds, col_greens, col_purples])[n]).append(col_tuple)

