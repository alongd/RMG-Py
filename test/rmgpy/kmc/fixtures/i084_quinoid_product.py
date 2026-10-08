# Exact stored product graph from the second intra_H_migration compile crash.

PRODUCT_GRAPH = """multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {7,S} {25,S}
2  C u0 p0 c0 {1,S} {4,S} {26,S} {27,S}
3  C u0 p0 c0 {1,S} {6,S} {28,S} {29,S}
4  C u0 p0 c0 {2,S} {8,S} {30,S} {31,S}
5  C u0 p0 c0 {6,S} {32,S} {33,S} {34,S}
6  C u0 p0 c0 {3,S} {5,S} {9,D}
7  C u0 p0 c0 {1,S} {12,S} {13,D}
8  C u0 p0 c0 {4,S} {10,S} {11,D}
9  C u0 p0 c0 {6,D} {14,S} {15,S}
10 C u0 p0 c0 {8,S} {16,D} {35,S}
11 C u0 p0 c0 {8,D} {18,S} {39,S}
12 C u0 p0 c0 {7,S} {19,D} {40,S}
13 C u0 p0 c0 {7,D} {21,S} {44,S}
14 C u0 p0 c0 {9,S} {22,D} {45,S}
15 C u0 p0 c0 {9,S} {23,D} {46,S}
16 C u0 p0 c0 {10,D} {17,S} {36,S}
17 C u0 p0 c0 {16,S} {18,D} {37,S}
18 C u0 p0 c0 {11,S} {17,D} {38,S}
19 C u0 p0 c0 {12,D} {20,S} {41,S}
20 C u0 p0 c0 {19,S} {21,D} {42,S}
21 C u0 p0 c0 {13,S} {20,D} {43,S}
22 C u0 p0 c0 {14,D} {24,S} {47,S}
23 C u0 p0 c0 {15,D} {24,S} {48,S}
24 C u2 p0 c0 {22,S} {23,S}
25 H u0 p0 c0 {1,S}
26 H u0 p0 c0 {2,S}
27 H u0 p0 c0 {2,S}
28 H u0 p0 c0 {3,S}
29 H u0 p0 c0 {3,S}
30 H u0 p0 c0 {4,S}
31 H u0 p0 c0 {4,S}
32 H u0 p0 c0 {5,S}
33 H u0 p0 c0 {5,S}
34 H u0 p0 c0 {5,S}
35 H u0 p0 c0 {10,S}
36 H u0 p0 c0 {16,S}
37 H u0 p0 c0 {17,S}
38 H u0 p0 c0 {18,S}
39 H u0 p0 c0 {11,S}
40 H u0 p0 c0 {12,S}
41 H u0 p0 c0 {19,S}
42 H u0 p0 c0 {20,S}
43 H u0 p0 c0 {21,S}
44 H u0 p0 c0 {13,S}
45 H u0 p0 c0 {14,S}
46 H u0 p0 c0 {15,S}
47 H u0 p0 c0 {22,S}
48 H u0 p0 c0 {23,S}
"""

# Complete stored transition captured from the bounded compiler path.  Its atom
# order is the compiler's canonical explicit-H order, shared with atom_map and
# bond_ops below.
CAPTURED_PRODUCT_GRAPH = """multiplicity 3
1  H u0 p0 c0 {26,S}
2  H u0 p0 c0 {27,S}
3  H u0 p0 c0 {28,S}
4  H u0 p0 c0 {29,S}
5  H u0 p0 c0 {32,S}
6  H u0 p0 c0 {33,S}
7  H u0 p0 c0 {34,S}
8  H u0 p0 c0 {35,S}
9  H u0 p0 c0 {36,S}
10 H u0 p0 c0 {37,S}
11 H u0 p0 c0 {38,S}
12 H u0 p0 c0 {39,S}
13 H u0 p0 c0 {40,S}
14 H u0 p0 c0 {41,S}
15 H u0 p0 c0 {44,S}
16 H u0 p0 c0 {44,S}
17 H u0 p0 c0 {44,S}
18 H u0 p0 c0 {45,S}
19 H u0 p0 c0 {45,S}
20 H u0 p0 c0 {46,S}
21 H u0 p0 c0 {46,S}
22 H u0 p0 c0 {47,S}
23 H u0 p0 c0 {47,S}
24 H u0 p0 c0 {48,S}
25 C u2 p0 c0 {26,S} {27,S}
26 C u0 p0 c0 {1,S} {25,S} {28,D}
27 C u0 p0 c0 {2,S} {25,S} {29,D}
28 C u0 p0 c0 {3,S} {26,D} {30,S}
29 C u0 p0 c0 {4,S} {27,D} {30,S}
30 C u0 p0 c0 {28,S} {29,S} {31,D}
31 C u0 p0 c0 {30,D} {44,S} {45,S}
32 C u0 p0 c0 {5,S} {34,B} {35,B}
33 C u0 p0 c0 {6,S} {36,B} {37,B}
34 C u0 p0 c0 {7,S} {32,B} {38,B}
35 C u0 p0 c0 {8,S} {32,B} {39,B}
36 C u0 p0 c0 {9,S} {33,B} {40,B}
37 C u0 p0 c0 {10,S} {33,B} {41,B}
38 C u0 p0 c0 {11,S} {34,B} {42,B}
39 C u0 p0 c0 {12,S} {35,B} {42,B}
40 C u0 p0 c0 {13,S} {36,B} {43,B}
41 C u0 p0 c0 {14,S} {37,B} {43,B}
42 C u0 p0 c0 {38,B} {39,B} {46,S}
43 C u0 p0 c0 {40,B} {41,B} {48,S}
44 C u0 p0 c0 {15,S} {16,S} {17,S} {31,S}
45 C u0 p0 c0 {18,S} {19,S} {31,S} {48,S}
46 C u0 p0 c0 {20,S} {21,S} {42,S} {47,S}
47 C u0 p0 c0 {22,S} {23,S} {46,S} {48,S}
48 C u0 p0 c0 {24,S} {43,S} {45,S} {47,S}
"""

ATOM_MAP = {
    0: 6,
    1: 23,
    2: 0,
    3: 8,
    4: 7,
    5: 1,
    6: 2,
    7: 11,
    8: 12,
    9: 9,
    10: 10,
    11: 3,
    12: 4,
    13: 15,
    14: 16,
    15: 13,
    16: 14,
    17: 5,
    18: 18,
    19: 17,
    20: 19,
    21: 20,
    22: 22,
    23: 21,
}

BOND_OPS = [
    {"action": "form", "atoms": [25, 0], "order": "1.0"},
    {"action": "break", "atoms": [24, 41], "order": "1.0"},
    {"action": "form", "atoms": [24, 41], "order": "2.0"},
    {"action": "break", "atoms": [41, 35], "order": "1.5"},
    {"action": "form", "atoms": [41, 35], "order": "1.0"},
    {"action": "break", "atoms": [41, 36], "order": "1.5"},
    {"action": "form", "atoms": [41, 36], "order": "1.0"},
    {"action": "break", "atoms": [35, 29], "order": "1.5"},
    {"action": "form", "atoms": [35, 29], "order": "2.0"},
    {"action": "break", "atoms": [36, 30], "order": "1.5"},
    {"action": "form", "atoms": [36, 30], "order": "2.0"},
    {"action": "break", "atoms": [29, 26], "order": "1.5"},
    {"action": "form", "atoms": [29, 26], "order": "1.0"},
    {"action": "break", "atoms": [26, 30], "order": "1.5"},
    {"action": "form", "atoms": [26, 30], "order": "1.0"},
    {"action": "break", "atoms": [26, 0], "order": "1.0"},
    {"action": "set_radical", "atom": 24, "value": 0},
    {"action": "set_radical", "atom": 25, "value": 0},
    {"action": "set_radical", "atom": 26, "value": 2},
]
