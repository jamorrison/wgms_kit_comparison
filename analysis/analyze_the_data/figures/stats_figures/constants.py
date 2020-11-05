"""Module to contain any values that won't change."""
# Map sample to short name
SAMPLE_KIT = {
    'FtubeAkapaBC': 'Kapa',
    'FtubeAneb': 'NEB',
    'FtubeApbat': 'PBAT',
    'FtubeAswift': 'Swift',
    'FtubeAkapaBCrep2': 'Kapa',
    'FtubeAnebRep2': 'NEB',
    'FtubeAswiftRep2': 'Swift',
    'FtubeAneb10ng': 'Low NEB',
    'FtubeAswift10ng': 'Low Swift',
    'FtubeAneb10ngRep2': 'Low NEB',
    'FtubeAswift10ngRep2': 'Low Swift',
    'FtubeBkapaBC': 'Kapa',
    'FtubeBneb': 'NEB',
    'FtubeBpbat': 'PBAT',
    'FtubeBswift': 'Swift',
    'FtubeBkapaBCrep2': 'Kapa',
    'FtubeBnebRep2': 'NEB',
    'FtubeBswiftRep2': 'Swift',
    'FtubeBneb10ng': 'Low NEB',
    'FtubeBswift10ng': 'Low Swift',
    'FtubeBneb10ngRep2': 'Low NEB',
    'FtubeBswift10ngRep2': 'Low Swift'
}

# Map sample to technical replicate
SAMPLE_REP = {
    'FtubeAkapaBC': '1',
    'FtubeAneb': '1',
    'FtubeApbat': '1',
    'FtubeAswift': '1',
    'FtubeAkapaBCrep2': '2',
    'FtubeAnebRep2': '2',
    'FtubeAswiftRep2': '2',
    'FtubeAneb10ng': '1',
    'FtubeAswift10ng': '1',
    'FtubeAneb10ngRep2': '2',
    'FtubeAswift10ngRep2': '2',
    'FtubeBkapaBC': '1',
    'FtubeBneb': '1',
    'FtubeBpbat': '1',
    'FtubeBswift': '1',
    'FtubeBkapaBCrep2': '2',
    'FtubeBnebRep2': '2',
    'FtubeBswiftRep2': '2',
    'FtubeBneb10ng': '1',
    'FtubeBswift10ng': '1',
    'FtubeBneb10ngRep2': '2',
    'FtubeBswift10ngRep2': '2'
}

# Map sample to short name + replicate
SAMPLE_NAMES = {
    'FtubeAkapaBC': 'Kapa Rep. 1',
    'FtubeAneb': 'NEB Rep. 1',
    'FtubeApbat': 'PBAT Rep. 1',
    'FtubeAswift': 'Swift Rep. 1',
    'FtubeAkapaBCrep2': 'Kapa Rep. 2',
    'FtubeAnebRep2': 'NEB Rep. 2',
    'FtubeAswiftRep2': 'Swift Rep. 2',
    'FtubeAneb10ng': 'Low NEB Rep. 1',
    'FtubeAswift10ng': 'Low Swift Rep. 1',
    'FtubeAneb10ngRep2': 'Low NEB Rep. 2',
    'FtubeAswift10ngRep2': 'Low Swift Rep. 2',
    'FtubeBkapaBC': 'Kapa Rep. 1',
    'FtubeBneb': 'NEB Rep. 1',
    'FtubeBpbat': 'PBAT Rep. 1',
    'FtubeBswift': 'Swift Rep. 1',
    'FtubeBkapaBCrep2': 'Kapa Rep. 2',
    'FtubeBnebRep2': 'NEB Rep. 2',
    'FtubeBswiftRep2': 'Swift Rep. 2',
    'FtubeBneb10ng': 'Low NEB Rep. 1',
    'FtubeBswift10ng': 'Low Swift Rep. 1',
    'FtubeBneb10ngRep2': 'Low NEB Rep. 2',
    'FtubeBswift10ngRep2': 'Low Swift Rep. 2'
}

# Map sample to biological sample
SAMPLE_GROUP = {
    'FtubeAkapaBC': 'A',
    'FtubeAneb': 'A',
    'FtubeApbat': 'A',
    'FtubeAswift': 'A',
    'FtubeAkapaBCrep2': 'A',
    'FtubeAnebRep2': 'A',
    'FtubeAswiftRep2': 'A',
    'FtubeAneb10ng': 'A',
    'FtubeAswift10ng': 'A',
    'FtubeAneb10ngRep2': 'A',
    'FtubeAswift10ngRep2': 'A',
    'FtubeBkapaBC': 'B',
    'FtubeBneb': 'B',
    'FtubeBpbat': 'B',
    'FtubeBswift': 'B',
    'FtubeBkapaBCrep2': 'B',
    'FtubeBnebRep2': 'B',
    'FtubeBswiftRep2': 'B',
    'FtubeBneb10ng': 'B',
    'FtubeBswift10ng': 'B',
    'FtubeBneb10ngRep2': 'B',
    'FtubeBswift10ngRep2': 'B'
}

# Map sample to what order to put them in on plots (biological sample order)
INTERLEAVE_ORDER = {
    'FtubeAkapaBC': 1,
    'FtubeBkapaBC': 2,
    'FtubeAkapaBCrep2': 3,
    'FtubeBkapaBCrep2': 4,
    'FtubeAneb': 5,
    'FtubeBneb': 6,
    'FtubeAnebRep2': 7,
    'FtubeBnebRep2': 8,
    'FtubeApbat': 9,
    'FtubeBpbat': 10,
    'FtubeAswift': 11,
    'FtubeBswift': 12,
    'FtubeAswiftRep2': 13,
    'FtubeBswiftRep2': 14,
    'FtubeAneb10ng': 15,
    'FtubeBneb10ng': 16,
    'FtubeAneb10ngRep2': 17,
    'FtubeBneb10ngRep2': 18,
    'FtubeAswift10ng': 19,
    'FtubeBswift10ng': 20,
    'FtubeAswift10ngRep2': 21,
    'FtubeBswift10ngRep2': 22
}

# Map sample to what order to put them on plots (replicate order)
REPLICATE_ORDER = {
    'FtubeAkapaBC': 1,
    'FtubeAkapaBCrep2': 2,
    'FtubeBkapaBC': 3,
    'FtubeBkapaBCrep2': 4,
    'FtubeAneb': 5,
    'FtubeAnebRep2': 6,
    'FtubeBneb': 7,
    'FtubeBnebRep2': 8,
    'FtubeApbat': 9,
    'FtubeBpbat': 10,
    'FtubeAswift': 11,
    'FtubeAswiftRep2': 12,
    'FtubeBswift': 13,
    'FtubeBswiftRep2': 14,
    'FtubeAneb10ng': 15,
    'FtubeAneb10ngRep2': 16,
    'FtubeBneb10ng': 17,
    'FtubeBneb10ngRep2': 18,
    'FtubeAswift10ng': 19,
    'FtubeAswift10ngRep2': 20,
    'FtubeBswift10ng': 21,
    'FtubeBswift10ngRep2': 22
}

# Map sample A samples to order on plots
A_ORDER = {
    'FtubeAkapaBC': 1,
    'FtubeAneb': 3,
    'FtubeApbat': 5,
    'FtubeAswift': 6,
    'FtubeAkapaBCrep2': 2,
    'FtubeAnebRep2': 4,
    'FtubeAswiftRep2': 7,
    'FtubeAneb10ng': 8,
    'FtubeAswift10ng': 10,
    'FtubeAneb10ngRep2': 9,
    'FtubeAswift10ngRep2': 11
}

# Map sample B samples to order on plots
B_ORDER = {
    'FtubeBkapaBC': 1,
    'FtubeBneb': 3,
    'FtubeBpbat': 5,
    'FtubeBswift': 6,
    'FtubeBkapaBCrep2': 2,
    'FtubeBnebRep2': 4,
    'FtubeBswiftRep2': 7,
    'FtubeBneb10ng': 8,
    'FtubeBswift10ng': 10,
    'FtubeBneb10ngRep2': 9,
    'FtubeBswift10ngRep2': 11
}

# Map sample to color on plot
SAMPLE_COLOR = {
    'FtubeAkapaBC': '#D81B60',
    'FtubeAneb': '#1E88E5',
    'FtubeApbat': '#FFC107',
    'FtubeAswift': '#004D40',
    'FtubeAkapaBCrep2': '#D81B60',
    'FtubeAnebRep2': '#1E88E5',
    'FtubeAswiftRep2': '#004D40',
    'FtubeAneb10ng': '#1E88E5',
    'FtubeAswift10ng': '#004D40',
    'FtubeAneb10ngRep2': '#1E88E5',
    'FtubeAswift10ngRep2': '#004D40',
    'FtubeBkapaBC': '#D81B60',
    'FtubeBneb': '#1E88E5',
    'FtubeBpbat': '#FFC107',
    'FtubeBswift': '#004D40',
    'FtubeBkapaBCrep2': '#D81B60',
    'FtubeBnebRep2': '#1E88E5',
    'FtubeBswiftRep2': '#004D40',
    'FtubeBneb10ng': '#1E88E5',
    'FtubeBswift10ng': '#004D40',
    'FtubeBneb10ngRep2': '#1E88E5',
    'FtubeBswift10ngRep2': '#004D40'
}

# Map sample to order on plot (replicate average plots)
SAMPLE_PLOT_VALUE = {
    'FtubeAkapaBC': 11,
    'FtubeAkapaBCrep2': 11,
    'FtubeBkapaBC': 10,
    'FtubeBkapaBCrep2': 10,
    'FtubeAneb': 9,
    'FtubeAnebRep2': 9,
    'FtubeBneb': 8,
    'FtubeBnebRep2': 8,
    'FtubeApbat': 7,
    'FtubeBpbat': 6,
    'FtubeAswift': 5,
    'FtubeAswiftRep2': 5,
    'FtubeBswift': 4,
    'FtubeBswiftRep2': 4,
    'FtubeAneb10ng': 3,
    'FtubeAneb10ngRep2': 3,
    'FtubeBneb10ng': 2,
    'FtubeBneb10ngRep2': 2,
    'FtubeAswift10ng': 1,
    'FtubeAswift10ngRep2': 1,
    'FtubeBswift10ng': 0,
    'FtubeBswift10ngRep2': 0
}

# Map sample to style to use on line graphs
SAMPLE_STYLE = {
    'FtubeAkapaBC': 'o-',
    'FtubeAneb': 'o-',
    'FtubeApbat': 'o-',
    'FtubeAswift': 'o-',
    'FtubeAkapaBCrep2': 's-',
    'FtubeAnebRep2': 's-',
    'FtubeAswiftRep2': 's-',
    'FtubeAneb10ng': 'o--',
    'FtubeAswift10ng': 'o--',
    'FtubeAneb10ngRep2': 's--',
    'FtubeAswift10ngRep2': 's--',
    'FtubeBkapaBC': 'o-',
    'FtubeBneb': 'o-',
    'FtubeBpbat': 'o-',
    'FtubeBswift': 'o-',
    'FtubeBkapaBCrep2': 's-',
    'FtubeBnebRep2': 's-',
    'FtubeBswiftRep2': 's-',
    'FtubeBneb10ng': 'o--',
    'FtubeBswift10ng': 'o--',
    'FtubeBneb10ngRep2': 's--',
    'FtubeBswift10ngRep2': 's--'
}

# Map sample to style to use on replicate average plots
REP_FORMAT = {
    '1': 'x',
    '2': 'o'
}
