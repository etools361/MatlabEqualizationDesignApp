function [DataAll, C0] = eq_data()
C0 = [2,5,8,10];
wTabNxC2 = {
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
];
[
3.564,3.807,3.937,3.985,3.998,4.000,4.000,4.000,4.000;
];
[
0.642,0.599,0.571,0.552,0.540,0.531,0.524,0.519,0.515;
];
[
2.268,2.610,3.024,3.437,3.751,3.920,3.982,3.997,4.000;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.914,0.836,0.783,0.748,0.725,0.710,0.700,0.692,0.687;
];
[
1.509,1.781,2.135,2.568,3.049,3.489,3.792,3.939,3.988;
1.779,1.971,2.248,2.623,3.068,3.494,3.793,3.939,3.988;
];
[
0.337,0.312,0.295,0.283,0.274,0.267,0.263,0.260,0.258;
1.088,0.981,0.908,0.858,0.825,0.803,0.789,0.779,0.773;
];
[
1.174,1.372,1.638,1.981,2.408,2.898,3.373,3.727,3.914;
1.581,1.698,1.871,2.126,2.481,2.926,3.381,3.728,3.914;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.552,0.509,0.479,0.457,0.442,0.430,0.422,0.416,0.412;
1.217,1.087,0.996,0.933,0.891,0.863,0.845,0.833,0.825;
];
[
0.939,1.104,1.322,1.604,1.960,2.398,2.898,3.380,3.734;
0.988,1.135,1.338,1.609,1.961,2.399,2.898,3.380,3.734;
1.482,1.560,1.677,1.852,2.110,2.473,2.926,3.387,3.736;
];
[
0.229,0.212,0.201,0.192,0.185,0.180,0.176,0.174,0.172;
0.707,0.647,0.606,0.577,0.556,0.541,0.529,0.521,0.516;
1.321,1.171,1.065,0.991,0.940,0.906,0.884,0.869,0.860;
];
[
0.796,0.932,1.113,1.348,1.647,2.025,2.486,2.997,3.466;
0.875,0.987,1.145,1.362,1.652,2.026,2.486,2.997,3.466;
1.426,1.482,1.566,1.692,1.881,2.161,2.549,3.019,3.471;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.398,0.368,0.346,0.331,0.319,0.310,0.304,0.299,0.295;
0.826,0.751,0.700,0.664,0.639,0.621,0.607,0.597,0.590;
1.408,1.241,1.121,1.037,0.979,0.940,0.914,0.897,0.886;
];
[
0.682,0.801,0.958,1.161,1.420,1.748,2.158,2.646,3.161;
0.700,0.812,0.964,1.163,1.421,1.748,2.158,2.646,3.161;
0.801,0.889,1.015,1.191,1.432,1.751,2.158,2.646,3.161;
1.390,1.433,1.496,1.591,1.735,1.952,2.269,2.693,3.175;
];
[
0.174,0.161,0.152,0.145,0.140,0.136,0.133,0.131,0.129;
0.530,0.487,0.457,0.436,0.421,0.409,0.400,0.393,0.388;
0.922,0.835,0.774,0.732,0.703,0.682,0.666,0.655,0.647;
1.483,1.301,1.170,1.076,1.012,0.967,0.938,0.918,0.906;
];
[
0.603,0.706,0.843,1.021,1.248,1.537,1.901,2.349,2.863;
0.632,0.725,0.854,1.025,1.250,1.537,1.901,2.349,2.863;
0.751,0.821,0.923,1.068,1.271,1.544,1.902,2.349,2.863;
1.367,1.400,1.450,1.524,1.637,1.808,2.064,2.430,2.893;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.312,0.288,0.271,0.259,0.250,0.243,0.237,0.233,0.230;
0.636,0.582,0.545,0.519,0.500,0.486,0.475,0.467,0.460;
1.004,0.904,0.834,0.786,0.753,0.730,0.712,0.700,0.691;
1.549,1.353,1.212,1.110,1.039,0.991,0.958,0.936,0.922;
];
[
0.536,0.628,0.751,0.910,1.113,1.371,1.696,2.104,2.593;
0.544,0.634,0.754,0.911,1.114,1.371,1.697,2.104,2.593;
0.582,0.661,0.770,0.919,1.116,1.371,1.697,2.104,2.593;
0.715,0.773,0.857,0.977,1.149,1.385,1.701,2.105,2.593;
1.351,1.377,1.417,1.477,1.568,1.706,1.915,2.225,2.645;
];
[
0.140,0.130,0.122,0.117,0.113,0.110,0.107,0.105,0.104;
0.425,0.391,0.368,0.351,0.338,0.329,0.321,0.316,0.311;
0.725,0.661,0.617,0.586,0.564,0.548,0.535,0.526,0.519;
1.074,0.963,0.885,0.832,0.794,0.769,0.750,0.736,0.726;
1.607,1.400,1.249,1.141,1.064,1.011,0.975,0.951,0.935;
];
[
0.485,0.568,0.679,0.822,1.005,1.237,1.531,1.902,2.359;
0.499,0.577,0.683,0.824,1.005,1.237,1.531,1.902,2.359;
0.545,0.612,0.706,0.836,1.010,1.239,1.532,1.902,2.359;
0.689,0.737,0.807,0.909,1.055,1.260,1.539,1.904,2.359;
1.339,1.360,1.393,1.442,1.517,1.631,1.805,2.066,2.438;
];
};
wTabNxC5 = {
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
];
[
8.749,9.416,9.802,9.954,9.993,9.999,10.000,10.000,10.000;
];
[
0.647,0.603,0.573,0.554,0.540,0.531,0.524,0.519,0.515;
];
[
5.597,6.401,7.384,8.411,9.254,9.749,9.943,9.991,9.999;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.921,0.841,0.787,0.751,0.727,0.711,0.700,0.692,0.687;
];
[
3.744,4.404,5.251,6.282,7.434,8.541,9.371,9.809,9.961;
4.432,4.891,5.545,6.426,7.488,8.555,9.373,9.809,9.961;
];
[
0.338,0.313,0.296,0.284,0.274,0.268,0.263,0.260,0.258;
1.096,0.987,0.912,0.861,0.827,0.805,0.790,0.780,0.773;
];
[
2.923,3.409,4.054,4.881,5.899,7.066,8.241,9.186,9.732;
3.955,4.236,4.651,5.254,6.090,7.142,8.262,9.189,9.732;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.554,0.510,0.480,0.458,0.442,0.431,0.423,0.417,0.412;
1.225,1.094,1.001,0.937,0.894,0.865,0.846,0.833,0.825;
];
[
2.341,2.749,3.285,3.971,4.828,5.873,7.065,8.257,9.205;
2.464,2.827,3.325,3.985,4.832,5.873,7.065,8.257,9.205;
3.716,3.906,4.187,4.604,5.215,6.068,7.141,8.278,9.209;
];
[
0.230,0.213,0.201,0.192,0.185,0.180,0.177,0.174,0.172;
0.708,0.648,0.606,0.577,0.557,0.541,0.530,0.522,0.516;
1.329,1.178,1.070,0.994,0.943,0.908,0.885,0.870,0.860;
];
[
1.986,2.323,2.769,3.346,4.076,4.984,6.079,7.304,8.478;
2.184,2.462,2.852,3.384,4.088,4.987,6.080,7.304,8.479;
3.579,3.716,3.919,4.222,4.674,5.335,6.247,7.364,8.493;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.399,0.368,0.346,0.331,0.319,0.311,0.304,0.299,0.295;
0.827,0.752,0.701,0.665,0.640,0.621,0.608,0.598,0.591;
1.416,1.248,1.127,1.041,0.982,0.942,0.915,0.898,0.886;
];
[
1.703,1.998,2.387,2.889,3.523,4.318,5.299,6.460,7.704;
1.747,2.026,2.401,2.893,3.524,4.318,5.299,6.460,7.705;
2.002,2.220,2.530,2.964,3.553,4.327,5.300,6.460,7.704;
3.493,3.596,3.750,3.980,4.326,4.841,5.587,6.585,7.744;
];
[
0.174,0.161,0.152,0.145,0.140,0.136,0.133,0.131,0.129;
0.530,0.487,0.457,0.436,0.421,0.409,0.400,0.393,0.388;
0.924,0.836,0.775,0.733,0.703,0.682,0.667,0.655,0.647;
1.491,1.307,1.175,1.081,1.015,0.970,0.939,0.919,0.906;
];
[
1.505,1.761,2.102,2.542,3.102,3.806,4.684,5.753,6.977;
1.578,1.810,2.129,2.553,3.105,3.807,4.684,5.753,6.977;
1.877,2.052,2.304,2.661,3.158,3.826,4.688,5.754,6.977;
3.436,3.517,3.637,3.818,4.090,4.499,5.106,5.965,7.059;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.312,0.288,0.271,0.259,0.250,0.243,0.238,0.233,0.230;
0.636,0.583,0.545,0.519,0.500,0.486,0.475,0.467,0.461;
1.005,0.905,0.835,0.787,0.754,0.730,0.713,0.700,0.691;
1.557,1.360,1.217,1.115,1.043,0.993,0.960,0.937,0.923;
];
[
1.338,1.569,1.874,2.269,2.770,3.402,4.192,5.170,6.332;
1.359,1.582,1.881,2.271,2.771,3.402,4.192,5.170,6.332;
1.455,1.650,1.922,2.291,2.777,3.403,4.193,5.170,6.332;
1.789,1.931,2.139,2.437,2.859,3.438,4.203,5.172,6.332;
3.396,3.461,3.558,3.703,3.923,4.255,4.754,5.482,6.471;
];
[
0.140,0.130,0.122,0.117,0.113,0.110,0.107,0.105,0.104;
0.425,0.391,0.368,0.351,0.339,0.329,0.322,0.316,0.311;
0.725,0.661,0.617,0.587,0.565,0.548,0.536,0.526,0.519;
1.075,0.964,0.886,0.832,0.795,0.769,0.751,0.737,0.727;
1.616,1.407,1.254,1.145,1.067,1.013,0.977,0.952,0.936;
];
[
1.212,1.419,1.694,2.049,2.502,3.074,3.793,4.687,5.775;
1.248,1.442,1.706,2.054,2.504,3.074,3.793,4.687,5.775;
1.362,1.528,1.762,2.084,2.516,3.078,3.793,4.687,5.774;
1.723,1.842,2.015,2.268,2.629,3.133,3.813,4.692,5.775;
3.367,3.420,3.500,3.620,3.801,4.075,4.490,5.108,5.984;
];
};
wTabNxC8 = {
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
];
[
13.376,14.562,15.425,15.851,15.976,15.998,16.000,16.000,16.000;
];
[
0.658,0.610,0.578,0.556,0.542,0.531,0.524,0.519,0.515;
];
[
8.744,9.880,11.275,12.805,14.236,15.289,15.817,15.972,15.997;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.935,0.852,0.795,0.756,0.730,0.713,0.701,0.693,0.687;
];
[
5.912,6.906,8.156,9.647,11.310,12.994,14.456,15.440,15.874;
7.052,7.721,8.657,9.900,11.410,13.022,14.462,15.441,15.874;
];
[
0.339,0.314,0.297,0.285,0.275,0.269,0.264,0.260,0.258;
1.111,0.999,0.921,0.868,0.832,0.808,0.792,0.781,0.773;
];
[
4.643,5.392,6.374,7.606,9.086,10.767,12.514,14.098,15.244;
6.342,6.761,7.369,8.236,9.419,10.906,12.557,14.107,15.244;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.556,0.511,0.481,0.459,0.444,0.432,0.424,0.417,0.413;
1.242,1.106,1.011,0.944,0.899,0.869,0.849,0.835,0.826;
];
[
3.726,4.365,5.196,6.241,7.522,9.043,10.760,12.536,14.132;
3.926,4.492,5.260,6.265,7.528,9.044,10.760,12.536,14.132;
5.983,6.268,6.686,7.298,8.178,9.383,10.900,12.578,14.141;
];
[
0.230,0.213,0.201,0.192,0.186,0.181,0.177,0.174,0.172;
0.710,0.650,0.608,0.579,0.558,0.543,0.531,0.523,0.517;
1.346,1.191,1.080,1.002,0.949,0.913,0.889,0.873,0.862;
];
[
3.167,3.696,4.395,5.287,6.397,7.747,9.337,11.103,12.883;
3.486,3.922,4.528,5.349,6.418,7.752,9.337,11.103,12.883;
5.776,5.982,6.287,6.737,7.398,8.346,9.633,11.217,12.915;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.399,0.368,0.347,0.331,0.320,0.311,0.305,0.300,0.296;
0.830,0.754,0.702,0.666,0.641,0.623,0.609,0.599,0.592;
1.433,1.261,1.137,1.049,0.988,0.947,0.919,0.900,0.888;
];
[
2.717,3.183,3.797,4.579,5.557,6.759,8.206,9.881,11.691;
2.789,3.228,3.819,4.587,5.559,6.759,8.206,9.881,11.691;
3.200,3.542,4.028,4.701,5.607,6.773,8.208,9.881,11.691;
5.646,5.802,6.034,6.378,6.889,7.638,8.701,10.107,11.769;
];
[
0.174,0.161,0.152,0.145,0.140,0.137,0.134,0.131,0.130;
0.531,0.488,0.458,0.437,0.421,0.410,0.401,0.394,0.389;
0.927,0.838,0.776,0.734,0.704,0.683,0.668,0.657,0.648;
1.508,1.321,1.185,1.089,1.021,0.975,0.943,0.922,0.908;
];
[
2.404,2.810,3.347,4.038,4.908,5.986,7.301,8.862,10.624;
2.521,2.889,3.391,4.056,4.913,5.987,7.301,8.862,10.624;
3.003,3.278,3.673,4.232,5.000,6.019,7.309,8.863,10.625;
5.559,5.682,5.864,6.135,6.540,7.141,8.017,9.232,10.775;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.312,0.288,0.272,0.259,0.250,0.243,0.238,0.234,0.231;
0.637,0.583,0.546,0.520,0.501,0.487,0.476,0.468,0.462;
1.008,0.907,0.837,0.788,0.755,0.731,0.714,0.702,0.692;
1.574,1.373,1.227,1.123,1.049,0.998,0.963,0.940,0.925;
];
[
2.137,2.504,2.988,3.610,4.393,5.369,6.568,8.014,9.694;
2.172,2.525,2.998,3.613,4.394,5.369,6.568,8.014,9.694;
2.325,2.635,3.065,3.645,4.405,5.371,6.568,8.014,9.694;
2.863,3.088,3.414,3.882,4.538,5.429,6.585,8.017,9.695;
5.499,5.597,5.744,5.963,6.292,6.784,7.510,8.550,9.942;
];
[
0.140,0.130,0.122,0.117,0.113,0.110,0.107,0.105,0.104;
0.425,0.391,0.368,0.351,0.339,0.329,0.322,0.316,0.312;
0.725,0.662,0.618,0.587,0.565,0.549,0.537,0.527,0.520;
1.077,0.966,0.888,0.834,0.796,0.770,0.752,0.738,0.728;
1.632,1.420,1.265,1.153,1.074,1.018,0.981,0.955,0.938;
];
[
1.934,2.266,2.701,3.264,3.975,4.864,5.964,7.304,8.890;
1.992,2.303,2.721,3.271,3.977,4.865,5.964,7.303,8.891;
2.176,2.441,2.812,3.320,3.997,4.870,5.965,7.303,8.891;
2.759,2.947,3.220,3.616,4.180,4.960,5.997,7.311,8.891;
5.454,5.535,5.656,5.837,6.110,6.518,7.128,8.019,9.255;
];
};
wTabNxC10 = {
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
];
[
15.814,17.180,18.348,19.187,19.677,19.901,19.978,19.997,20.000;
];
[
0.670,0.619,0.584,0.560,0.544,0.533,0.525,0.519,0.515;
];
[
10.700,11.951,13.468,15.135,16.759,18.126,19.091,19.645,19.894;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.950,0.863,0.803,0.762,0.735,0.716,0.703,0.694,0.687;
];
[
7.301,8.477,9.926,11.611,13.453,15.316,17.012,18.356,19.249;
8.778,9.541,10.590,11.957,13.597,15.362,17.023,18.358,19.249;
];
[
0.341,0.316,0.298,0.286,0.276,0.270,0.265,0.261,0.258;
1.127,1.012,0.931,0.875,0.837,0.812,0.795,0.783,0.774;
];
[
5.764,6.670,7.840,9.282,10.972,12.844,14.774,16.575,18.050;
7.953,8.439,9.136,10.114,11.422,13.041,14.840,16.592,18.053;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.558,0.513,0.482,0.461,0.445,0.433,0.425,0.418,0.414;
1.258,1.119,1.021,0.952,0.905,0.873,0.852,0.837,0.827;
];
[
4.637,5.418,6.426,7.677,9.179,10.917,12.831,14.793,16.613;
4.889,5.579,6.508,7.708,9.187,10.918,12.831,14.794,16.614;
7.530,7.864,8.351,9.053,10.045,11.378,13.030,14.859,16.629;
];
[
0.231,0.213,0.201,0.193,0.186,0.181,0.178,0.175,0.173;
0.713,0.652,0.609,0.580,0.559,0.544,0.533,0.524,0.518;
1.363,1.204,1.090,1.010,0.955,0.918,0.892,0.875,0.864;
];
[
3.946,4.598,5.452,6.533,7.859,9.436,11.244,13.206,15.177;
4.349,4.883,5.622,6.612,7.886,9.442,11.244,13.206,15.178;
7.286,7.529,7.886,8.409,9.165,10.230,11.648,13.371,15.229;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.400,0.369,0.347,0.332,0.321,0.312,0.305,0.300,0.296;
0.832,0.756,0.704,0.668,0.642,0.624,0.610,0.600,0.593;
1.450,1.274,1.147,1.057,0.995,0.952,0.923,0.903,0.890;
];
[
3.388,3.965,4.719,5.675,6.856,8.283,9.959,11.849,13.852;
3.479,4.021,4.747,5.685,6.858,8.283,9.959,11.849,13.852;
3.996,4.417,5.011,5.831,6.920,8.301,9.963,11.850,13.852;
7.133,7.317,7.590,7.992,8.582,9.435,10.621,12.162,13.969;
];
[
0.175,0.161,0.152,0.146,0.141,0.137,0.134,0.132,0.130;
0.532,0.488,0.459,0.438,0.422,0.410,0.401,0.395,0.390;
0.929,0.840,0.778,0.735,0.706,0.685,0.669,0.658,0.649;
1.526,1.334,1.196,1.097,1.028,0.980,0.947,0.925,0.911;
];
[
2.999,3.502,4.166,5.015,6.074,7.368,8.915,10.704,12.672;
3.147,3.601,4.221,5.038,6.080,7.369,8.915,10.704,12.672;
3.753,4.092,4.577,5.259,6.191,7.410,8.925,10.705,12.672;
7.030,7.174,7.389,7.708,8.180,8.871,9.860,11.205,12.886;
];
[
0.000,0.000,-0.000,0.000,0.000,-0.000,0.000,0.000,-0.000;
0.313,0.289,0.272,0.260,0.251,0.244,0.238,0.234,0.231;
0.638,0.584,0.547,0.521,0.502,0.488,0.477,0.468,0.462;
1.011,0.909,0.839,0.790,0.756,0.732,0.715,0.703,0.693;
1.592,1.387,1.238,1.131,1.056,1.003,0.967,0.943,0.927;
];
[
2.668,3.123,3.722,4.488,5.447,6.628,8.056,9.736,11.636;
2.711,3.149,3.735,4.493,5.448,6.628,8.056,9.736,11.637;
2.903,3.287,3.819,4.533,5.462,6.631,8.056,9.736,11.637;
3.580,3.858,4.260,4.833,5.632,6.705,8.078,9.741,11.637;
6.957,7.074,7.248,7.506,7.891,8.460,9.289,10.452,11.980;
];
[
0.140,0.130,0.122,0.117,0.113,0.110,0.107,0.106,0.104;
0.425,0.392,0.368,0.352,0.339,0.330,0.322,0.317,0.312;
0.726,0.663,0.619,0.588,0.566,0.550,0.537,0.528,0.521;
1.080,0.969,0.890,0.835,0.797,0.771,0.753,0.739,0.729;
1.650,1.434,1.276,1.161,1.080,1.024,0.985,0.958,0.941;
];
[
2.415,2.827,3.367,4.062,4.936,6.019,7.340,8.915,10.732;
2.488,2.874,3.392,4.072,4.939,6.019,7.340,8.915,10.732;
2.719,3.047,3.506,4.133,4.964,6.026,7.341,8.915,10.732;
3.452,3.684,4.021,4.507,5.196,6.141,7.382,8.925,10.734;
6.905,7.001,7.144,7.358,7.677,8.153,8.854,9.859,11.229;
];
};
DataAll = {wTabNxC2, wTabNxC5, wTabNxC8, wTabNxC10};