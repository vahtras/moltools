import unittest, os
import numpy as np

import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr

from molecules import Water
from pdbreader import System, Residue, Atom
from property import Property

WATER_FILE_TARGZ = os.path.join(os.path.dirname(__file__), 'hfqua_water.tar.gz')
PRO_1_FILE_TARGZ = os.path.join(os.path.dirname(__file__), 'b3lyp_A-PRO1-ready.tar.gz')
PRO_2_FILE_TARGZ = os.path.join(os.path.dirname(__file__), 'b3lyp_A-PRO2-ready.tar.gz')
#PRO_3_FILE_TARGZ = os.path.join(os.path.dirname(__file__), 'b3lyp_A-GLY3-ready.tar.gz')
CON_1_FILE_TARGZ = os.path.join(os.path.dirname(__file__), 'b3lyp_A-PRO1-concap.tar.gz')
#CON_2_FILE_TARGZ = os.path.join(os.path.dirname(__file__), 'b3lyp_A-PRO2-concap.tar.gz')

_ppg_string = """ATOM  0     N    PRO A   1      -1.849  -2.684  38.880                       N
ATOM  1     H1   PRO A   1      -1.416  -1.763  38.665                       H
ATOM  2     CD   PRO A   1      -3.089  -2.462  39.690                       C
ATOM  3     HD1  PRO A   1      -2.996  -1.543  40.304                       H
ATOM  4     HD2  PRO A   1      -3.193  -3.314  40.393                       H
ATOM  5     CG   PRO A   1      -4.284  -2.400  38.715                       C
ATOM  6     HG1  PRO A   1      -4.562  -1.351  38.505                       H
ATOM  7     HG2  PRO A   1      -5.177  -2.882  39.142                       H
ATOM  8     CB   PRO A   1      -3.815  -3.104  37.431                       C
ATOM  9     HB1  PRO A   1      -4.022  -2.483  36.539                       H
ATOM  10    HB2  PRO A   1      -4.349  -4.055  37.278                       H
ATOM  11    CA   PRO A   1      -2.294  -3.354  37.617                       C
ATOM  12    HA   PRO A   1      -2.091  -4.447  37.753                       H
ATOM  13    C    PRO A   1      -1.516  -2.817  36.415                       C
ATOM  14    O    PRO A   1      -0.887  -1.757  36.460                       O
ATOM  0     N    PRO A   2      -1.476  -3.576  35.256                       N
ATOM  1     CD   PRO A   2      -2.189  -4.861  35.009                       C
ATOM  2     HD1  PRO A   2      -3.285  -4.689  34.939                       H
ATOM  3     HD2  PRO A   2      -2.001  -5.567  35.842                       H
ATOM  4     CG   PRO A   2      -1.617  -5.370  33.672                       C
ATOM  5     HG1  PRO A   2      -2.402  -5.845  33.060                       H
ATOM  6     HG2  PRO A   2      -0.842  -6.138  33.843                       H
ATOM  7     CB   PRO A   2      -1.016  -4.147  32.953                       C
ATOM  8     HB1  PRO A   2      -1.635  -3.864  32.079                       H
ATOM  9     HB2  PRO A   2      -0.007  -4.355  32.541                       H
ATOM  10    CA   PRO A   2      -0.952  -3.007  33.987                       C
ATOM  11    HA   PRO A   2       0.107  -2.647  34.163                       H
ATOM  12    C    PRO A   2      -1.824  -1.812  33.536                       C
ATOM  13    O    PRO A   2      -2.939  -1.565  33.980                       O
ATOM  0     N    GLY A   3      -1.204  -1.013  32.594                       N
ATOM  1     H    GLY A   3      -0.231  -1.176  32.307                       H
ATOM  2     CA   GLY A   3      -1.832   0.206  32.079                       C
ATOM  3     HA1  GLY A   3      -1.016   0.964  31.867                       H
ATOM  4     HA2  GLY A   3      -2.485   0.671  32.863                       H
ATOM  5     C    GLY A   3      -2.643  -0.034  30.812                       C
ATOM  6     O    GLY A   3      -2.917  -1.122  30.320                       O
ATOM  0     N    PRO A   4      -3.126   1.135  30.213                       N
ATOM  1     CD   PRO A   4      -2.771   2.520  30.648                       C
ATOM  2     HD1  PRO A   4      -1.682   2.603  30.875                       H
ATOM  3     HD2  PRO A   4      -3.314   2.754  31.590                       H
ATOM  4     CG   PRO A   4      -3.199   3.441  29.496                       C
ATOM  5     HG1  PRO A   4      -2.314   3.778  28.921                       H
ATOM  6     HG2  PRO A   4      -3.696   4.351  29.870                       H
ATOM  7     CB   PRO A   4      -4.137   2.611  28.601                       C
ATOM  8     HB1  PRO A   4      -3.976   2.857  27.532                       H
ATOM  9     HB2  PRO A   4      -5.194   2.835  28.826                       H
ATOM  10    CA   PRO A   4      -3.828   1.128  28.908                       C
ATOM  11    HA   PRO A   4      -4.755   0.507  29.010                       H
ATOM  12    C    PRO A   4      -2.943   0.550  27.783                       C
ATOM  13    O    PRO A   4      -1.729   0.385  27.881                       O
ATOM  0     N    PRO A   5      -3.587   0.187  26.608                       N
ATOM  1     CD   PRO A   5      -5.045   0.298  26.321                       C
ATOM  2     HD1  PRO A   5      -5.348   1.365  26.251                       H
ATOM  3     HD2  PRO A   5      -5.629  -0.177  27.136                       H
ATOM  4     CG   PRO A   5      -5.240  -0.424  24.974                       C
ATOM  5     HG1  PRO A   5      -5.982   0.097  24.347                       H
ATOM  6     HG2  PRO A   5      -5.623  -1.448  25.129                       H
ATOM  7     CB   PRO A   5      -3.861  -0.458  24.291                       C
ATOM  8     HB1  PRO A   5      -3.822   0.268  23.456                       H
ATOM  9     HB2  PRO A   5      -3.644  -1.447  23.831                       H
ATOM  10    CA   PRO A   5      -2.818  -0.122  25.374                       C
ATOM  11    HA   PRO A   5      -2.134  -0.999  25.591                       H
ATOM  12    C    PRO A   5      -1.966   1.095  24.947                       C
ATOM  13    O    PRO A   5      -2.119   2.235  25.368                       O
ATOM  0     N    GLY A   6      -0.962   0.764  24.053                       N
ATOM  1     H    GLY A   6      -0.798  -0.206  23.755                       H
ATOM  2     CA   GLY A   6      -0.022   1.763  23.543                       C
ATOM  3     HA1  GLY A   6       0.972   1.250  23.360                       H
ATOM  4     HA2  GLY A   6       0.169   2.550  24.317                       H
ATOM  5     C    GLY A   6      -0.504   2.432  22.259                       C
ATOM  6     O    GLY A   6      -1.630   2.359  21.782                       O
ATOM  0     N    PRO A   7       0.460   3.232  21.635                       N
ATOM  1     CD   PRO A   7       1.896   3.301  22.040                       C
ATOM  2     HD1  PRO A   7       2.298   2.284  22.255                       H
ATOM  3     HD2  PRO A   7       1.979   3.891  22.979                       H
ATOM  4     CG   PRO A   7       2.627   3.981  20.873                       C
ATOM  5     HG1  PRO A   7       3.198   3.234  20.286                       H
ATOM  6     HG2  PRO A   7       3.357   4.724  21.231                       H
ATOM  7     CB   PRO A   7       1.540   4.633  19.999                       C
ATOM  8     HB1  PRO A   7       1.800   4.548  18.925                       H
ATOM  9     HB2  PRO A   7       1.450   5.709  20.223                       H
ATOM  10    CA   PRO A   7       0.218   3.903  20.336                       C
ATOM  11    HA   PRO A   7      -0.643   4.610  20.459                       H
ATOM  12    C    PRO A   7      -0.101   2.894  19.213                       C
ATOM  13    O    PRO A   7       0.091   1.683  19.302                       O
ATOM  0     N    PRO A   8      -0.655   3.412  18.049                       N
ATOM  1     CD   PRO A   8      -0.976   4.841  17.776                       C
ATOM  2     HD1  PRO A   8      -0.045   5.443  17.704                       H
ATOM  3     HD2  PRO A   8      -1.595   5.254  18.598                       H
ATOM  4     CG   PRO A   8      -1.732   4.827  16.434                       C
ATOM  5     HG1  PRO A   8      -1.455   5.694  15.811                       H
ATOM  6     HG2  PRO A   8      -2.823   4.893  16.596                       H
ATOM  7     CB   PRO A   8      -1.367   3.504  15.738                       C
ATOM  8     HB1  PRO A   8      -0.668   3.685  14.899                       H
ATOM  9     HB2  PRO A   8      -2.253   3.011  15.280                       H
ATOM  10    CA   PRO A   8      -0.733   2.596  16.809                       C
ATOM  11    HA   PRO A   8      -1.369   1.683  17.022                       H
ATOM  12    C    PRO A   8       0.679   2.143  16.370                       C
ATOM  13    O    PRO A   8       1.725   2.621  16.791                       O
ATOM  0     N    GLY A   9       0.652   1.095  15.466                       N
ATOM  1     H    GLY A   9      -0.228   0.657  15.167                       H
ATOM  2     CA   GLY A   9       1.881   0.500  14.939                       C
ATOM  3     HA1  GLY A   9       1.685  -0.600  14.748                       H
ATOM  4     HA2  GLY A   9       2.696   0.544  15.707                       H
ATOM  5     C    GLY A   9       2.365   1.171  13.658                       C
ATOM  6     O    GLY A   9       1.957   2.230  13.197                       O
ATOM  0     N    PRO A  10       3.409   0.496  13.016                       N
ATOM  1     CD   PRO A  10       3.904  -0.860  13.400                       C
ATOM  2     HD1  PRO A  10       3.053  -1.549  13.613                       H
ATOM  3     HD2  PRO A  10       4.498  -0.777  14.336                       H
ATOM  4     CG   PRO A  10       4.760  -1.342  12.220                       C
ATOM  5     HG1  PRO A  10       4.213  -2.102  11.629                       H
ATOM  6     HG2  PRO A  10       5.691  -1.822  12.565                       H
ATOM  7     CB   PRO A  10       5.052  -0.101  11.359                       C
ATOM  8     HB1  PRO A  10       5.039  -0.362  10.281                       H
ATOM  9     HB2  PRO A  10       6.055   0.303  11.578                       H
ATOM  10    CA   PRO A  10       3.966   0.941  11.716                       C
ATOM  11    HA   PRO A  10       4.385   1.972  11.848                       H
ATOM  12    C    PRO A  10       2.897   0.958  10.603                       C
ATOM  13    O    PRO A  10       1.799   0.414  10.697                       O
ATOM  0     N    PRO A  11       3.217   1.652   9.444                       N
ATOM  1     CD   PRO A  11       4.484   2.386   9.167                       C
ATOM  2     HD1  PRO A  11       5.335   1.676   9.081                       H
ATOM  3     HD2  PRO A  11       4.701   3.092   9.993                       H
ATOM  4     CG   PRO A  11       4.235   3.116   7.834                       C
ATOM  5     HG1  PRO A  11       5.140   3.115   7.203                       H
ATOM  6     HG2  PRO A  11       3.975   4.176   8.008                       H
ATOM  7     CB   PRO A  11       3.075   2.381   7.140                       C
ATOM  8     HB1  PRO A  11       3.449   1.774   6.292                       H
ATOM  9     HB2  PRO A  11       2.337   3.084   6.695                       H
ATOM  10    CA   PRO A  11       2.405   1.496   8.208                       C
ATOM  11    HA   PRO A  11       1.346   1.830   8.433                       H
ATOM  12    C    PRO A  11       2.388   0.017   7.756                       C
ATOM  13    O    PRO A  11       3.159  -0.843   8.164                       O
ATOM  0     N    GLY A  12       1.373  -0.261   6.857                       N
ATOM  1     H    GLY A  12       0.691   0.452   6.570                       H
ATOM  2     CA   GLY A  12       1.167  -1.607   6.318                       C
ATOM  3     HA1  GLY A  12       0.057  -1.745   6.134                       H
ATOM  4     HA2  GLY A  12       1.457  -2.378   7.076                       H
ATOM  5     C    GLY A  12       1.943  -1.855   5.029                       C
ATOM  6     O    GLY A  12       2.829  -1.144   4.570                       O
ATOM  0     N    PRO A  13       1.606  -3.047   4.378                       N
ATOM  1     CD   PRO A  13       0.463  -3.928   4.762                       C
ATOM  2     HD1  PRO A  13      -0.448  -3.324   4.985                       H
ATOM  3     HD2  PRO A  13       0.726  -4.478   5.691                       H
ATOM  4     CG   PRO A  13       0.251  -4.879   3.574                       C
ATOM  5     HG1  PRO A  13      -0.642  -4.579   2.992                       H
ATOM  6     HG2  PRO A  13       0.074  -5.913   3.911                       H
ATOM  7     CB   PRO A  13       1.518  -4.778   2.706                       C
ATOM  8     HB1  PRO A  13       1.258  -4.835   1.630                       H
ATOM  9     HB2  PRO A  13       2.205  -5.616   2.914                       H
ATOM  10    CA   PRO A  13       2.188  -3.433   3.071                       C
ATOM  11    HA   PRO A  13       3.299  -3.525   3.194                       H
ATOM  12    C    PRO A  13       1.877  -2.398   1.969                       C
ATOM  13    O    PRO A  13       1.030  -1.514   2.078                       O
ATOM  0     N    PRO A  14       2.626  -2.485   0.803                       N
ATOM  1     CD   PRO A  14       3.703  -3.472   0.509                       C
ATOM  2     HD1  PRO A  14       3.280  -4.495   0.417                       H
ATOM  3     HD2  PRO A  14       4.448  -3.475   1.331                       H
ATOM  4     CG   PRO A  14       4.316  -3.003  -0.824                       C
ATOM  5     HG1  PRO A  14       4.582  -3.861  -1.464                       H
ATOM  6     HG2  PRO A  14       5.251  -2.440  -0.651                       H
ATOM  7     CB   PRO A  14       3.262  -2.110  -1.502                       C
ATOM  8     HB1  PRO A  14       2.789  -2.641  -2.350                       H
ATOM  9     HB2  PRO A  14       3.709  -1.192  -1.943                       H
ATOM  10    CA   PRO A  14       2.226  -1.745  -0.422                       C
ATOM  11    HA   PRO A  14       2.230  -0.636  -0.187                       H
ATOM  12    C    PRO A  14       0.806  -2.165  -0.868                       C
ATOM  13    O    PRO A  14       0.217  -3.161  -0.465                       O
ATOM  0     N    GLY A  15       0.232  -1.271  -1.755                       N
ATOM  1     H    GLY A  15       0.707  -0.405  -2.037                       H
ATOM  2     CA   GLY A  15      -1.119  -1.469  -2.284                       C
ATOM  3     HA1  GLY A  15      -1.581  -0.449  -2.457                       H
ATOM  4     HA2  GLY A  15      -1.764  -1.981  -1.525                       H
ATOM  5     C    GLY A  15      -1.136  -2.274  -3.580                       C
ATOM  6     O    GLY A  15      -0.199  -2.903  -4.054                       O
ATOM  0     N    PRO A  16      -2.381  -2.302  -4.220                       N
ATOM  1     CD   PRO A  16      -3.559  -1.477  -3.818                       C
ATOM  2     HD1  PRO A  16      -3.253  -0.429  -3.589                       H
ATOM  3     HD2  PRO A  16      -3.997  -1.900  -2.888                       H
ATOM  4     CG   PRO A  16      -4.542  -1.550  -4.996                       C
ATOM  5     HG1  PRO A  16      -4.527  -0.604  -5.572                       H
ATOM  6     HG2  PRO A  16      -5.578  -1.691  -4.650                       H
ATOM  7     CB   PRO A  16      -4.076  -2.722  -5.878                       C
ATOM  8     HB1  PRO A  16      -4.219  -2.484  -6.951                       H
ATOM  9     HB2  PRO A  16      -4.668  -3.630  -5.669                       H
ATOM  10    CA   PRO A  16      -2.588  -2.962  -5.530                       C
ATOM  11    HA   PRO A  16      -2.342  -4.050  -5.419                       H
ATOM  12    C    PRO A  16      -1.705  -2.346  -6.635                       C
ATOM  13    O    PRO A  16      -1.112  -1.275  -6.523                       O
ATOM  0     N    PRO A  17      -1.575  -3.077  -7.809                       N
ATOM  1     CD   PRO A  17      -2.197  -4.397  -8.108                       C
ATOM  2     HD1  PRO A  17      -3.301  -4.299  -8.192                       H
ATOM  3     HD2  PRO A  17      -1.972  -5.115  -7.294                       H
ATOM  4     CG   PRO A  17      -1.576  -4.832  -9.449                       C
ATOM  5     HG1  PRO A  17      -2.320  -5.338 -10.088                       H
ATOM  6     HG2  PRO A  17      -0.758  -5.556  -9.288                       H
ATOM  7     CB   PRO A  17      -1.044  -3.555 -10.122                       C
ATOM  8     HB1  PRO A  17      -1.700  -3.254 -10.962                       H
ATOM  9     HB2  PRO A  17      -0.039  -3.703 -10.573                       H
ATOM  10    CA   PRO A  17      -0.995  -2.465  -9.033                       C
ATOM  11    HA   PRO A  17       0.066  -2.143  -8.803                       H
ATOM  12    C    PRO A  17      -1.820  -1.230  -9.464                       C
ATOM  13    O    PRO A  17      -2.942  -0.964  -9.051                       O
ATOM  0     N    GLY A  18      -1.141  -0.410 -10.349                       N
ATOM  1     H    GLY A  18      -0.175  -0.605 -10.638                       H
ATOM  2     CA   GLY A  18      -1.730   0.829 -10.860                       C
ATOM  3     HA1  GLY A  18      -0.892   1.571 -11.034                       H
ATOM  4     HA2  GLY A  18      -2.400   1.289 -10.088                       H
ATOM  5     C    GLY A  18      -2.519   0.625 -12.149                       C
ATOM  6     O    GLY A  18      -2.849  -0.450 -12.635                       O
ATOM  0     N    PRO A  19      -2.918   1.814 -12.770                       N
ATOM  1     CD   PRO A  19      -2.473   3.178 -12.354                       C
ATOM  2     HD1  PRO A  19      -1.380   3.192 -12.135                       H
ATOM  3     HD2  PRO A  19      -2.998   3.460 -11.415                       H
ATOM  4     CG   PRO A  19      -2.844   4.111 -13.518                       C
ATOM  5     HG1  PRO A  19      -1.942   4.381 -14.100                       H
ATOM  6     HG2  PRO A  19      -3.280   5.056 -13.155                       H
ATOM  7     CB   PRO A  19      -3.838   3.332 -14.398                       C
ATOM  8     HB1  PRO A  19      -3.667   3.554 -15.470                       H
ATOM  9     HB2  PRO A  19      -4.877   3.627 -14.171                       H
ATOM  10    CA   PRO A  19      -3.624   1.835 -14.073                       C
ATOM  11    HA   PRO A  19      -4.590   1.277 -13.959                       H
ATOM  12    C    PRO A  19      -2.786   1.188 -15.195                       C
ATOM  13    O    PRO A  19      -1.587   0.937 -15.099                       O
ATOM  0     N    PRO A  20      -3.458   0.863 -16.366                       N
ATOM  1     CD   PRO A  20      -4.907   1.070 -16.648                       C
ATOM  2     HD1  PRO A  20      -5.139   2.155 -16.721                       H
ATOM  3     HD2  PRO A  20      -5.518   0.638 -15.830                       H
ATOM  4     CG   PRO A  20      -5.154   0.358 -17.991                       C
ATOM  5     HG1  PRO A  20      -5.865   0.923 -18.617                       H
ATOM  6     HG2  PRO A  20      -5.602  -0.639 -17.831                       H
ATOM  7     CB   PRO A  20      -3.784   0.233 -18.680                       C
ATOM  8     HB1  PRO A  20      -3.699   0.958 -19.513                       H
ATOM  9     HB2  PRO A  20      -3.635  -0.766 -19.144                       H
ATOM  10    CA   PRO A  20      -2.717   0.496 -17.601                       C
ATOM  11    HA   PRO A  20      -2.098  -0.429 -17.383                       H
ATOM  12    C    PRO A  20      -1.777   1.646 -18.031                       C
ATOM  13    O    PRO A  20      -1.843   2.794 -17.610                       O
ATOM  0     N    GLY A  21      -0.807   1.239 -18.931                       N
ATOM  1     H    GLY A  21      -0.715   0.259 -19.225                       H
ATOM  2     CA   GLY A  21       0.214   2.161 -19.433                       C
ATOM  3     HA1  GLY A  21       1.161   1.568 -19.622                       H
ATOM  4     HA2  GLY A  21       0.473   2.920 -18.650                       H
ATOM  5     C    GLY A  21      -0.207   2.882 -20.708                       C
ATOM  6     O    GLY A  21      -1.331   2.900 -21.195                       O
ATOM  0     N    PRO A  22       0.819   3.613 -21.317                       N
ATOM  1     CD   PRO A  22       2.253   3.571 -20.901                       C
ATOM  2     HD1  PRO A  22       2.578   2.524 -20.694                       H
ATOM  3     HD2  PRO A  22       2.373   4.142 -19.955                       H
ATOM  4     CG   PRO A  22       3.041   4.207 -22.056                       C
ATOM  5     HG1  PRO A  22       3.552   3.426 -22.653                       H
ATOM  6     HG2  PRO A  22       3.827   4.884 -21.685                       H
ATOM  7     CB   PRO A  22       2.012   4.957 -22.920                       C
ATOM  8     HB1  PRO A  22       2.271   4.879 -23.994                       H
ATOM  9     HB2  PRO A  22       1.999   6.032 -22.671                       H
ATOM  10    CA   PRO A  22       0.638   4.320 -22.607                       C
ATOM  11    HA   PRO A  22      -0.172   5.086 -22.478                       H
ATOM  12    C    PRO A  22       0.253   3.353 -23.746                       C
ATOM  13    O    PRO A  22       0.354   2.130 -23.674                       O
ATOM  0     N    PRO A  23      -0.253   3.927 -24.905                       N
ATOM  1     CD   PRO A  23      -0.468   5.378 -25.162                       C
ATOM  2     HD1  PRO A  23       0.506   5.908 -25.240                       H
ATOM  3     HD2  PRO A  23      -1.045   5.831 -24.330                       H
ATOM  4     CG   PRO A  23      -1.238   5.434 -26.495                       C
ATOM  5     HG1  PRO A  23      -0.906   6.286 -27.112                       H
ATOM  6     HG2  PRO A  23      -2.319   5.577 -26.319                       H
ATOM  7     CB   PRO A  23      -0.977   4.096 -27.209                       C
ATOM  8     HB1  PRO A  23      -0.264   4.232 -28.044                       H
ATOM  9     HB2  PRO A  23      -1.896   3.678 -27.673                       H
ATOM  10    CA   PRO A  23      -0.417   3.129 -26.148                       C
ATOM  11    HA   PRO A  23      -1.131   2.277 -25.932                       H
ATOM  12    C    PRO A  23       0.942   2.548 -26.603                       C
ATOM  13    O    PRO A  23       2.034   2.925 -26.200                       O
ATOM  0     N    GLY A  24       0.803   1.507 -27.507                       N
ATOM  1     H    GLY A  24      -0.120   1.150 -27.783                       H
ATOM  2     CA   GLY A  24       1.960   0.768 -28.015                       C
ATOM  3     HA1  GLY A  24       1.646  -0.307 -28.186                       H
ATOM  4     HA2  GLY A  24       2.773   0.738 -27.244                       H
ATOM  5     C    GLY A  24       2.518   1.349 -29.309                       C
ATOM  6     O    GLY A  24       2.225   2.428 -29.808                       O
ATOM  0     N    PRO A  25       3.484   0.543 -29.923                       N
ATOM  1     CD   PRO A  25       3.841  -0.840 -29.483                       C
ATOM  2     HD1  PRO A  25       2.927  -1.422 -29.215                       H
ATOM  3     HD2  PRO A  25       4.466  -0.778 -28.565                       H
ATOM  4     CG   PRO A  25       4.609  -1.464 -30.657                       C
ATOM  5     HG1  PRO A  25       3.959  -2.166 -31.215                       H
ATOM  6     HG2  PRO A  25       5.476  -2.049 -30.307                       H
ATOM  7     CB   PRO A  25       5.043  -0.295 -31.558                       C
ATOM  8     HB1  PRO A  25       4.996  -0.582 -32.626                       H
ATOM  9     HB2  PRO A  25       6.090  -0.010 -31.350                       H
ATOM  10    CA   PRO A  25       4.096   0.883 -31.228                       C
ATOM  11    HA   PRO A  25       4.646   1.854 -31.115                       H
ATOM  12    C    PRO A  25       3.038   1.024 -32.341                       C
ATOM  13    O    PRO A  25       1.888   0.593 -32.265                       O
ATOM  0     N    PRO A  26       3.416   1.720 -33.481                       N
ATOM  1     CD   PRO A  26       4.782   2.225 -33.802                       C
ATOM  2     HD1  PRO A  26       5.548   1.443 -33.626                       H
ATOM  3     HD2  PRO A  26       5.011   3.090 -33.144                       H
ATOM  4     CG   PRO A  26       4.741   2.637 -35.286                       C
ATOM  5     HG1  PRO A  26       5.330   1.930 -35.900                       H
ATOM  6     HG2  PRO A  26       5.188   3.633 -35.440                       H
ATOM  7     CB   PRO A  26       3.264   2.613 -35.717                       C
ATOM  8     HB1  PRO A  26       3.146   2.196 -36.739                       H
ATOM  9     HB2  PRO A  26       2.839   3.631 -35.753                       H
ATOM  10    CA   PRO A  26       2.520   1.767 -34.665                       C
ATOM  11    HA   PRO A  26       1.539   2.243 -34.362                       H
ATOM  12    C    PRO A  26       2.267   0.353 -35.230                       C
ATOM  13    O    PRO A  26       2.981  -0.619 -35.035                       O
ATOM  0     N    GLY A  27       1.083   0.261 -35.966                       N
ATOM  1     H    GLY A  27       0.539   1.109 -36.185                       H
ATOM  2     CA   GLY A  27       0.896  -0.856 -36.901                       C
ATOM  3     HA1  GLY A  27      -0.207  -1.006 -37.077                       H
ATOM  4     HA2  GLY A  27       1.268  -1.806 -36.440                       H
ATOM  5     C    GLY A  27       1.610  -0.565 -38.218                       C
ATOM  6     O    GLY A  27       2.145   0.509 -38.468                       O
ATOM  0     N    PRO A  28       1.556  -1.560 -39.205                       N
ATOM  1     CD   PRO A  28       1.239  -2.993 -38.914                       C
ATOM  2     HD1  PRO A  28       0.161  -3.162 -39.127                       H
ATOM  3     HD2  PRO A  28       1.415  -3.235 -37.845                       H
ATOM  4     CG   PRO A  28       2.137  -3.834 -39.842                       C
ATOM  5     HG1  PRO A  28       1.520  -4.371 -40.589                       H
ATOM  6     HG2  PRO A  28       2.700  -4.597 -39.281                       H
ATOM  7     CB   PRO A  28       3.081  -2.849 -40.556                       C
ATOM  8     HB1  PRO A  28       3.217  -3.130 -41.617                       H
ATOM  9     HB2  PRO A  28       4.082  -2.847 -40.092                       H
ATOM  10    CA   PRO A  28       2.422  -1.460 -40.416                       C
ATOM  11    HA   PRO A  28       3.186  -0.651 -40.288                       H
ATOM  12    C    PRO A  28       1.515  -1.183 -41.636                       C
ATOM  13    O    PRO A  28       0.949  -2.086 -42.249                       O
ATOM  0     N    PRO A  29       1.409   0.115 -42.102                       N
ATOM  1     CD   PRO A  29       1.897   1.355 -41.431                       C
ATOM  2     HD1  PRO A  29       2.964   1.513 -41.692                       H
ATOM  3     HD2  PRO A  29       1.826   1.266 -40.320                       H
ATOM  4     CG   PRO A  29       1.011   2.493 -41.961                       C
ATOM  5     HG1  PRO A  29       1.612   3.381 -42.217                       H
ATOM  6     HG2  PRO A  29       0.288   2.820 -41.183                       H
ATOM  7     CB   PRO A  29       0.263   1.952 -43.191                       C
ATOM  8     HB1  PRO A  29       0.594   2.451 -44.123                       H
ATOM  9     HB2  PRO A  29      -0.826   2.146 -43.109                       H
ATOM  10    CA   PRO A  29       0.527   0.438 -43.262                       C
ATOM  11    HA   PRO A  29      -0.439  -0.139 -43.187                       H
ATOM  12    C    PRO A  29       1.219  -0.003 -44.561                       C
ATOM  13    OC1  PRO A  29       1.113   0.624 -45.589                       O
ATOM  14    HCO  PRO A  29       1.793  -0.941 -44.493                       H
ATOM  0     N    PRO B   1      -2.862   1.060  41.428                       N
ATOM  1     H1   PRO B   1      -1.974   0.586  41.179                       H
ATOM  2     CD   PRO B   1      -2.563   2.198  42.357                       C
ATOM  3     HD1  PRO B   1      -1.614   2.012  42.906                       H
ATOM  4     HD2  PRO B   1      -3.381   2.249  43.103                       H
ATOM  5     CG   PRO B   1      -2.491   3.489  41.517                       C
ATOM  6     HG1  PRO B   1      -1.438   3.734  41.281                       H
ATOM  7     HG2  PRO B   1      -2.900   4.352  42.062                       H
ATOM  8     CB   PRO B   1      -3.283   3.193  40.232                       C
ATOM  9     HB1  PRO B   1      -2.740   3.557  39.340                       H
ATOM  10    HB2  PRO B   1      -4.256   3.709  40.238                       H
ATOM  11    CA   PRO B   1      -3.486   1.655  40.202                       C
ATOM  12    HA   PRO B   1      -4.575   1.399  40.251                       H
ATOM  13    C    PRO B   1      -2.859   1.068  38.933                       C
ATOM  14    O    PRO B   1      -1.760   0.512  38.944                       O
ATOM  0     N    PRO B   2      -3.577   1.145  37.751                       N
ATOM  1     CD   PRO B   2      -4.878   1.844  37.552                       C
ATOM  2     HD1  PRO B   2      -4.737   2.945  37.624                       H
ATOM  3     HD2  PRO B   2      -5.601   1.535  38.333                       H
ATOM  4     CG   PRO B   2      -5.335   1.430  36.141                       C
ATOM  5     HG1  PRO B   2      -5.817   2.272  35.616                       H
ATOM  6     HG2  PRO B   2      -6.082   0.618  36.192                       H
ATOM  7     CB   PRO B   2      -4.076   0.956  35.392                       C
ATOM  8     HB1  PRO B   2      -3.765   1.703  34.637                       H
ATOM  9     HB2  PRO B   2      -4.255   0.013  34.832                       H
ATOM  10    CA   PRO B   2      -2.977   0.743  36.452                       C
ATOM  11    HA   PRO B   2      -2.669  -0.340  36.522                       H
ATOM  12    C    PRO B   2      -1.746   1.624  36.136                       C
ATOM  13    O    PRO B   2      -1.612   2.782  36.505                       O
ATOM  0     N    GLY B   3      -0.781   0.965  35.390                       N
ATOM  1     H    GLY B   3      -0.849  -0.038  35.176                       H
ATOM  2     CA   GLY B   3       0.470   1.617  34.996                       C
ATOM  3     HA1  GLY B   3       1.272   0.822  34.920                       H
ATOM  4     HA2  GLY B   3       0.808   2.332  35.802                       H
ATOM  5     C    GLY B   3       0.351   2.372  33.678                       C
ATOM  6     O    GLY B   3      -0.690   2.704  33.124                       O
ATOM  0     N    PRO B   4       1.577   2.746  33.114                       N
ATOM  1     CD   PRO B   4       2.912   2.280  33.597                       C
ATOM  2     HD1  PRO B   4       2.885   1.196  33.863                       H
ATOM  3     HD2  PRO B   4       3.176   2.834  34.523                       H
ATOM  4     CG   PRO B   4       3.893   2.577  32.453                       C
ATOM  5     HG1  PRO B   4       4.136   1.649  31.900                       H
ATOM  6     HG2  PRO B   4       4.847   2.977  32.833                       H
ATOM  7     CB   PRO B   4       3.187   3.585  31.528                       C
ATOM  8     HB1  PRO B   4       3.436   3.386  30.467                       H
ATOM  9     HB2  PRO B   4       3.515   4.615  31.746                       H
ATOM  10    CA   PRO B   4       1.672   3.437  31.806                       C
ATOM  11    HA   PRO B   4       1.154   4.427  31.895                       H
ATOM  12    C    PRO B   4       1.028   2.620  30.666                       C
ATOM  13    O    PRO B   4       0.760   1.423  30.743                       O
ATOM  0     N    PRO B   5       0.732   3.308  29.498                       N
ATOM  1     CD   PRO B   5       0.971   4.754  29.230                       C
ATOM  2     HD1  PRO B   5       2.061   4.960  29.158                       H
ATOM  3     HD2  PRO B   5       0.553   5.368  30.054                       H
ATOM  4     CG   PRO B   5       0.264   5.031  27.889                       C
ATOM  5     HG1  PRO B   5       0.848   5.730  27.268                       H
ATOM  6     HG2  PRO B   5      -0.720   5.503  28.054                       H
ATOM  7     CB   PRO B   5       0.101   3.670  27.190                       C
ATOM  8     HB1  PRO B   5       0.818   3.574  26.351                       H
ATOM  9     HB2  PRO B   5      -0.905   3.550  26.732                       H
ATOM  10    CA   PRO B   5       0.343   2.586  28.258                       C
ATOM  11    HA   PRO B   5      -0.595   1.987  28.474                       H
ATOM  12    C    PRO B   5       1.469   1.626  27.814                       C
ATOM  13    O    PRO B   5       2.623   1.667  28.224                       O
ATOM  0     N    GLY B   6       1.034   0.663  26.920                       N
ATOM  1     H    GLY B   6       0.051   0.594  26.630                       H
ATOM  2     CA   GLY B   6       1.935  -0.365  26.394                       C
ATOM  3     HA1  GLY B   6       1.327  -1.304  26.209                       H
ATOM  4     HA2  GLY B   6       2.706  -0.636  27.160                       H
ATOM  5     C    GLY B   6       2.637   0.060  25.109                       C
ATOM  6     O    GLY B   6       2.672   1.192  24.643                       O
ATOM  0     N    PRO B   7       3.333  -0.972  24.468                       N
ATOM  1     CD   PRO B   7       3.262  -2.412  24.860                       C
ATOM  2     HD1  PRO B   7       2.212  -2.713  25.082                       H
ATOM  3     HD2  PRO B   7       3.849  -2.563  25.792                       H
ATOM  4     CG   PRO B   7       3.854  -3.195  23.678                       C
ATOM  5     HG1  PRO B   7       3.048  -3.681  23.093                       H
ATOM  6     HG2  PRO B   7       4.523  -4.001  24.022                       H
ATOM  7     CB   PRO B   7       4.605  -2.170  22.810                       C
ATOM  8     HB1  PRO B   7       4.486  -2.408  21.734                       H
ATOM  9     HB2  PRO B   7       5.687  -2.191  23.024                       H
ATOM  10    CA   PRO B   7       4.016  -0.785  23.166                       C
ATOM  11    HA   PRO B   7       4.806  -0.001  23.293                       H
ATOM  12    C    PRO B   7       3.036  -0.352  22.055                       C
ATOM  13    O    PRO B   7       1.813  -0.425  22.151                       O
ATOM  0     N    PRO B   8       3.599   0.162  20.894                       N
ATOM  1     CD   PRO B   8       5.051   0.343  20.615                       C
ATOM  2     HD1  PRO B   8       5.558  -0.642  20.526                       H
ATOM  3     HD2  PRO B   8       5.529   0.908  21.441                       H
ATOM  4     CG   PRO B   8       5.105   1.116  19.283                       C
ATOM  5     HG1  PRO B   8       5.936   0.762  18.650                       H
ATOM  6     HG2  PRO B   8       5.280   2.192  19.459                       H
ATOM  7     CB   PRO B   8       3.748   0.894  18.592                       C
ATOM  8     HB1  PRO B   8       3.853   0.191  17.743                       H
ATOM  9     HB2  PRO B   8       3.342   1.830  18.150                       H
ATOM  10    CA   PRO B   8       2.788   0.338  19.662                       C
ATOM  11    HA   PRO B   8       1.944   1.059  19.890                       H
ATOM  12    C    PRO B   8       2.193  -1.015  19.208                       C
ATOM  13    O    PRO B   8       2.568  -2.109  19.611                       O
ATOM  0     N    GLY B   9       1.147  -0.872  18.313                       N
ATOM  1     H    GLY B   9       0.796   0.051  18.030                       H
ATOM  2     CA   GLY B   9       0.427  -2.028  17.777                       C
ATOM  3     HA1  GLY B   9      -0.649  -1.720  17.598                       H
ATOM  4     HA2  GLY B   9       0.395  -2.853  18.534                       H
ATOM  5     C    GLY B   9       1.037  -2.560  16.484                       C
ATOM  6     O    GLY B   9       2.128  -2.255  16.019                       O
ATOM  0     N    PRO B  10       0.254  -3.523  15.835                       N
ATOM  1     CD   PRO B  10      -1.141  -3.883  16.226                       C
ATOM  2     HD1  PRO B  10      -1.739  -2.970  16.456                       H
ATOM  3     HD2  PRO B  10      -1.111  -4.496  17.153                       H
ATOM  4     CG   PRO B  10      -1.717  -4.671  15.040                       C
ATOM  5     HG1  PRO B  10      -2.422  -4.041  14.462                       H
ATOM  6     HG2  PRO B  10      -2.286  -5.553  15.377                       H
ATOM  7     CB   PRO B  10      -0.518  -5.075  14.164                       C
ATOM  8     HB1  PRO B  10      -0.784  -5.022  13.089                       H
ATOM  9     HB2  PRO B  10      -0.216  -6.117  14.367                       H
ATOM  10    CA   PRO B  10       0.631  -4.105  14.525                       C
ATOM  11    HA   PRO B  10       1.615  -4.628  14.643                       H
ATOM  12    C    PRO B  10       0.748  -3.029  13.425                       C
ATOM  13    O    PRO B  10       0.318  -1.883  13.537                       O
ATOM  0     N    PRO B  11       1.398  -3.402  12.256                       N
ATOM  1     CD   PRO B  11       1.998  -4.732  11.958                       C
ATOM  2     HD1  PRO B  11       1.206  -5.507  11.868                       H
ATOM  3     HD2  PRO B  11       2.685  -5.030  12.776                       H
ATOM  4     CG   PRO B  11       2.741  -4.541  10.622                       C
ATOM  5     HG1  PRO B  11       2.644  -5.433   9.981                       H
ATOM  6     HG2  PRO B  11       3.822  -4.392  10.790                       H
ATOM  7     CB   PRO B  11       2.122  -3.304   9.948                       C
ATOM  8     HB1  PRO B  11       1.475  -3.604   9.102                       H
ATOM  9     HB2  PRO B  11       2.892  -2.635   9.507                       H
ATOM  10    CA   PRO B  11       1.316  -2.562  11.032                       C
ATOM  11    HA   PRO B  11       1.756  -1.545  11.267                       H
ATOM  12    C    PRO B  11      -0.157  -2.391  10.592                       C
ATOM  13    O    PRO B  11      -1.087  -3.077  10.998                       O
ATOM  0     N    GLY B  12      -0.337  -1.342   9.709                       N
ATOM  1     H    GLY B  12       0.438  -0.731   9.424                       H
ATOM  2     CA   GLY B  12      -1.659  -0.994   9.184                       C
ATOM  3     HA1  GLY B  12      -1.686   0.126   9.014                       H
ATOM  4     HA2  GLY B  12      -2.451  -1.215   9.945                       H
ATOM  5     C    GLY B  12      -1.994  -1.724   7.887                       C
ATOM  6     O    GLY B  12      -1.380  -2.670   7.411                       O
ATOM  0     N    PRO B  13      -3.150  -1.259   7.250                       N
ATOM  1     CD   PRO B  13      -3.908  -0.038   7.656                       C
ATOM  2     HD1  PRO B  13      -3.214   0.804   7.886                       H
ATOM  3     HD2  PRO B  13      -4.476  -0.257   8.586                       H
ATOM  4     CG   PRO B  13      -4.842   0.284   6.479                       C
ATOM  5     HG1  PRO B  13      -4.456   1.150   5.905                       H
ATOM  6     HG2  PRO B  13      -5.850   0.561   6.827                       H
ATOM  7     CB   PRO B  13      -4.875  -0.974   5.594                       C
ATOM  8     HB1  PRO B  13      -4.913  -0.696   4.521                       H
ATOM  9     HB2  PRO B  13      -5.777  -1.576   5.800                       H
ATOM  10    CA   PRO B  13      -3.603  -1.782   5.939                       C
ATOM  11    HA   PRO B  13      -3.806  -2.879   6.049                       H
ATOM  12    C    PRO B  13      -2.550  -1.563   4.833                       C
ATOM  13    O    PRO B  13      -1.583  -0.812   4.945                       O
ATOM  0     N    PRO B  14      -2.720  -2.284   3.659                       N
ATOM  1     CD   PRO B  14      -3.813  -3.252   3.360                       C
ATOM  2     HD1  PRO B  14      -4.789  -2.727   3.281                       H
ATOM  3     HD2  PRO B  14      -3.885  -4.004   4.172                       H
ATOM  4     CG   PRO B  14      -3.418  -3.891   2.015                       C
ATOM  5     HG1  PRO B  14      -4.303  -4.060   1.379                       H
ATOM  6     HG2  PRO B  14      -2.951  -4.880   2.171                       H
ATOM  7     CB   PRO B  14      -2.428  -2.924   1.343                       C
ATOM  8     HB1  PRO B  14      -2.914  -2.388   0.505                       H
ATOM  9     HB2  PRO B  14      -1.563  -3.455   0.889                       H
ATOM  10    CA   PRO B  14      -1.952  -1.944   2.432                       C
ATOM  11    HA   PRO B  14      -0.848  -2.063   2.657                       H
ATOM  12    C    PRO B  14      -2.231  -0.484   2.007                       C
ATOM  13    O    PRO B  14      -3.159   0.197   2.426                       O
ATOM  0     N    GLY B  15      -1.290   0.010   1.121                       N
ATOM  1     H    GLY B  15      -0.478  -0.545   0.825                       H
ATOM  2     CA   GLY B  15      -1.355   1.380   0.610                       C
ATOM  3     HA1  GLY B  15      -0.295   1.741   0.434                       H
ATOM  4     HA2  GLY B  15      -1.795   2.063   1.381                       H
ATOM  5     C    GLY B  15      -2.163   1.494  -0.679                       C
ATOM  6     O    GLY B  15      -2.887   0.631  -1.159                       O
ATOM  0     N    PRO B  16      -2.070   2.743  -1.305                       N
ATOM  1     CD   PRO B  16      -1.128   3.827  -0.896                       C
ATOM  2     HD1  PRO B  16      -0.115   3.414  -0.679                       H
ATOM  3     HD2  PRO B  16      -1.498   4.295   0.042                       H
ATOM  4     CG   PRO B  16      -1.110   4.827  -2.063                       C
ATOM  5     HG1  PRO B  16      -0.174   4.724  -2.646                       H
ATOM  6     HG2  PRO B  16      -1.144   5.868  -1.704                       H
ATOM  7     CB   PRO B  16      -2.329   4.492  -2.940                       C
ATOM  8     HB1  PRO B  16      -2.085   4.623  -4.014                       H
ATOM  9     HB2  PRO B  16      -3.171   5.170  -2.719                       H
ATOM  10    CA   PRO B  16      -2.716   3.032  -2.607                       C
ATOM  11    HA   PRO B  16      -3.822   2.896  -2.489                       H
ATOM  12    C    PRO B  16      -2.200   2.104  -3.727                       C
ATOM  13    O    PRO B  16      -1.194   1.406  -3.631                       O
ATOM  0     N    PRO B  17      -2.949   2.063  -4.896                       N
ATOM  1     CD   PRO B  17      -4.202   2.819  -5.177                       C
ATOM  2     HD1  PRO B  17      -3.993   3.908  -5.248                       H
ATOM  3     HD2  PRO B  17      -4.934   2.657  -4.360                       H
ATOM  4     CG   PRO B  17      -4.707   2.262  -6.522                       C
ATOM  5     HG1  PRO B  17      -5.139   3.062  -7.147                       H
ATOM  6     HG2  PRO B  17      -5.508   1.519  -6.365                       H
ATOM  7     CB   PRO B  17      -3.493   1.613  -7.211                       C
ATOM  8     HB1  PRO B  17      -3.135   2.246  -8.046                       H
ATOM  9     HB2  PRO B  17      -3.746   0.634  -7.673                       H
ATOM  10    CA   PRO B  17      -2.407   1.442  -6.133                       C
ATOM  11    HA   PRO B  17      -2.190   0.350  -5.918                       H
ATOM  12    C    PRO B  17      -1.099   2.145  -6.563                       C
ATOM  13    O    PRO B  17      -0.720   3.230  -6.138                       O
ATOM  0     N    GLY B  18      -0.357   1.400  -7.463                       N
ATOM  1     H    GLY B  18      -0.649   0.462  -7.763                       H
ATOM  2     CA   GLY B  18       0.930   1.872  -7.978                       C
ATOM  3     HA1  GLY B  18       1.585   0.967  -8.169                       H
ATOM  4     HA2  GLY B  18       1.459   2.485  -7.204                       H
ATOM  5     C    GLY B  18       0.793   2.691  -9.257                       C
ATOM  6     O    GLY B  18      -0.248   3.132  -9.727                       O
ATOM  0     N    PRO B  19       2.010   2.979  -9.886                       N
ATOM  1     CD   PRO B  19       3.327   2.397  -9.490                       C
ATOM  2     HD1  PRO B  19       3.234   1.305  -9.283                       H
ATOM  3     HD2  PRO B  19       3.668   2.879  -8.548                       H
ATOM  4     CG   PRO B  19       4.281   2.688 -10.658                       C
ATOM  5     HG1  PRO B  19       4.456   1.770 -11.253                       H
ATOM  6     HG2  PRO B  19       5.268   3.025 -10.301                       H
ATOM  7     CB   PRO B  19       3.596   3.764 -11.520                       C
ATOM  8     HB1  PRO B  19       3.789   3.583 -12.596                       H
ATOM  9     HB2  PRO B  19       3.994   4.766 -11.287                       H
ATOM  10    CA   PRO B  19       2.089   3.695 -11.181                       C
ATOM  11    HA   PRO B  19       1.630   4.710 -11.051                       H
ATOM  12    C    PRO B  19       1.350   2.939 -12.306                       C
ATOM  13    O    PRO B  19       0.982   1.770 -12.220                       O
ATOM  0     N    PRO B  20       1.083   3.654 -13.466                       N
ATOM  1     CD   PRO B  20       1.431   5.078 -13.734                       C
ATOM  2     HD1  PRO B  20       2.533   5.202 -13.814                       H
ATOM  3     HD2  PRO B  20       1.068   5.719 -12.906                       H
ATOM  4     CG   PRO B  20       0.737   5.410 -15.069                       C
ATOM  5     HG1  PRO B  20       1.366   6.068 -15.692                       H
ATOM  6     HG2  PRO B  20      -0.209   5.953 -14.896                       H
ATOM  7     CB   PRO B  20       0.470   4.067 -15.771                       C
ATOM  8     HB1  PRO B  20       1.175   3.921 -16.612                       H
ATOM  9     HB2  PRO B  20      -0.543   4.023 -16.225                       H
ATOM  10    CA   PRO B  20       0.636   2.966 -14.705                       C
ATOM  11    HA   PRO B  20      -0.343   2.438 -14.487                       H
ATOM  12    C    PRO B  20       1.687   1.926 -15.157                       C
ATOM  13    O    PRO B  20       2.838   1.875 -14.744                       O
ATOM  0     N    GLY B  21       1.182   1.009 -16.063                       N
ATOM  1     H    GLY B  21       0.196   1.017 -16.351                       H
ATOM  2     CA   GLY B  21       2.000  -0.084 -16.592                       C
ATOM  3     HA1  GLY B  21       1.321  -0.970 -16.787                       H
ATOM  4     HA2  GLY B  21       2.742  -0.422 -15.823                       H
ATOM  5     C    GLY B  21       2.741   0.292 -17.871                       C
ATOM  6     O    GLY B  21       2.858   1.418 -18.337                       O
ATOM  0     N    PRO B  22       3.365  -0.787 -18.508                       N
ATOM  1     CD   PRO B  22       3.194  -2.218 -18.117                       C
ATOM  2     HD1  PRO B  22       2.124  -2.448 -17.905                       H
ATOM  3     HD2  PRO B  22       3.760  -2.407 -17.179                       H
ATOM  4     CG   PRO B  22       3.745  -3.040 -19.292                       C
ATOM  5     HG1  PRO B  22       2.915  -3.470 -19.886                       H
ATOM  6     HG2  PRO B  22       4.355  -3.888 -18.941                       H
ATOM  7     CB   PRO B  22       4.573  -2.067 -20.151                       C
ATOM  8     HB1  PRO B  22       4.456  -2.300 -21.228                       H
ATOM  9     HB2  PRO B  22       5.648  -2.157 -19.919                       H
ATOM  10    CA   PRO B  22       4.069  -0.647 -19.805                       C
ATOM  11    HA   PRO B  22       4.906   0.087 -19.675                       H
ATOM  12    C    PRO B  22       3.125  -0.156 -20.923                       C
ATOM  13    O    PRO B  22       1.899  -0.146 -20.832                       O
ATOM  0     N    PRO B  23       3.725   0.313 -22.084                       N
ATOM  1     CD   PRO B  23       5.187   0.397 -22.360                       C
ATOM  2     HD1  PRO B  23       5.624  -0.620 -22.454                       H
ATOM  3     HD2  PRO B  23       5.700   0.921 -21.528                       H
ATOM  4     CG   PRO B  23       5.295   1.174 -23.685                       C
ATOM  5     HG1  PRO B  23       6.103   0.770 -24.319                       H
ATOM  6     HG2  PRO B  23       5.542   2.235 -23.501                       H
ATOM  7     CB   PRO B  23       3.928   1.048 -24.381                       C
ATOM  8     HB1  PRO B  23       3.987   0.341 -25.230                       H
ATOM  9     HB2  PRO B  23       3.591   2.011 -24.823                       H
ATOM  10    CA   PRO B  23       2.929   0.560 -23.315                       C
ATOM  11    HA   PRO B  23       2.144   1.343 -23.080                       H
ATOM  12    C    PRO B  23       2.229  -0.738 -23.777                       C
ATOM  13    O    PRO B  23       2.509  -1.862 -23.379                       O
ATOM  0     N    GLY B  24       1.203  -0.506 -24.678                       N
ATOM  1     H    GLY B  24       0.925   0.444 -24.952                       H
ATOM  2     CA   GLY B  24       0.378  -1.596 -25.202                       C
ATOM  3     HA1  GLY B  24      -0.663  -1.190 -25.393                       H
ATOM  4     HA2  GLY B  24       0.261  -2.402 -24.431                       H
ATOM  5     C    GLY B  24       0.937  -2.206 -26.483                       C
ATOM  6     O    GLY B  24       2.044  -1.997 -26.964                       O
ATOM  0     N    PRO B  25       0.075  -3.113 -27.110                       N
ATOM  1     CD   PRO B  25      -1.341  -3.358 -26.703                       C
ATOM  2     HD1  PRO B  25      -1.867  -2.398 -26.491                       H
ATOM  3     HD2  PRO B  25      -1.352  -3.950 -25.762                       H
ATOM  4     CG   PRO B  25      -1.982  -4.126 -27.869                       C
ATOM  5     HG1  PRO B  25      -2.639  -3.458 -28.460                       H
ATOM  6     HG2  PRO B  25      -2.616  -4.953 -27.509                       H
ATOM  7     CB   PRO B  25      -0.822  -4.643 -28.738                       C
ATOM  8     HB1  PRO B  25      -1.090  -4.604 -29.813                       H
ATOM  9     HB2  PRO B  25      -0.595  -5.698 -28.503                       H
ATOM  10    CA   PRO B  25       0.397  -3.749 -28.409                       C
ATOM  11    HA   PRO B  25       1.344  -4.336 -28.292                       H
ATOM  12    C    PRO B  25       0.572  -2.705 -29.530                       C
ATOM  13    O    PRO B  25       0.223  -1.529 -29.435                       O
ATOM  0     N    PRO B  26       1.176  -3.136 -30.705                       N
ATOM  1     CD   PRO B  26       1.670  -4.509 -31.006                       C
ATOM  2     HD1  PRO B  26       0.824  -5.229 -31.048                       H
ATOM  3     HD2  PRO B  26       2.369  -4.840 -30.211                       H
ATOM  4     CG   PRO B  26       2.367  -4.392 -32.375                       C
ATOM  5     HG1  PRO B  26       2.145  -5.268 -33.008                       H
ATOM  6     HG2  PRO B  26       3.465  -4.363 -32.257                       H
ATOM  7     CB   PRO B  26       1.858  -3.094 -33.029                       C
ATOM  8     HB1  PRO B  26       1.199  -3.320 -33.890                       H
ATOM  9     HB2  PRO B  26       2.685  -2.485 -33.450                       H
ATOM  10    CA   PRO B  26       1.108  -2.308 -31.937                       C
ATOM  11    HA   PRO B  26       1.595  -1.307 -31.722                       H
ATOM  12    C    PRO B  26      -0.366  -2.088 -32.345                       C
ATOM  13    O    PRO B  26      -1.311  -2.748 -31.929                       O
ATOM  0     N    GLY B  27      -0.550  -1.008 -33.194                       N
ATOM  1     H    GLY B  27       0.236  -0.461 -33.561                       H
ATOM  2     CA   GLY B  27      -1.883  -0.671 -33.701                       C
ATOM  3     HA1  GLY B  27      -1.990   0.451 -33.713                       H
ATOM  4     HA2  GLY B  27      -2.665  -1.053 -32.992                       H
ATOM  5     C    GLY B  27      -2.138  -1.272 -35.078                       C
ATOM  6     O    GLY B  27      -1.453  -2.123 -35.627                       O
ATOM  0     N    PRO B  28      -3.301  -0.820 -35.715                       N
ATOM  1     CD   PRO B  28      -4.183   0.283 -35.239                       C
ATOM  2     HD1  PRO B  28      -3.592   1.207 -35.045                       H
ATOM  3     HD2  PRO B  28      -4.661  -0.014 -34.280                       H
ATOM  4     CG   PRO B  28      -5.225   0.488 -36.353                       C
ATOM  5     HG1  PRO B  28      -5.008   1.423 -36.919                       H
ATOM  6     HG2  PRO B  28      -6.242   0.596 -35.944                       H
ATOM  7     CB   PRO B  28      -5.114  -0.736 -37.280                       C
ATOM  8     HB1  PRO B  28      -5.299  -0.465 -38.333                       H
ATOM  9     HB2  PRO B  28      -5.874  -1.495 -37.014                       H
ATOM  10    CA   PRO B  28      -3.697  -1.308 -37.057                       C
ATOM  11    HA   PRO B  28      -3.681  -2.431 -37.050                       H
ATOM  12    C    PRO B  28      -2.723  -0.799 -38.141                       C
ATOM  13    O    PRO B  28      -1.903   0.098 -37.955                       O
ATOM  0     N    PRO B  29      -2.832  -1.361 -39.403                       N
ATOM  1     CD   PRO B  29      -3.658  -2.551 -39.760                       C
ATOM  2     HD1  PRO B  29      -4.701  -2.223 -39.959                       H
ATOM  3     HD2  PRO B  29      -3.681  -3.282 -38.926                       H
ATOM  4     CG   PRO B  29      -3.009  -3.146 -41.021                       C
ATOM  5     HG1  PRO B  29      -3.771  -3.451 -41.759                       H
ATOM  6     HG2  PRO B  29      -2.432  -4.056 -40.775                       H
ATOM  7     CB   PRO B  29      -2.085  -2.061 -41.596                       C
ATOM  8     HB1  PRO B  29      -2.520  -1.607 -42.510                       H
ATOM  9     HB2  PRO B  29      -1.104  -2.483 -41.910                       H
ATOM  10    CA   PRO B  29      -1.892  -0.995 -40.505                       C
ATOM  11    HA   PRO B  29      -0.833  -0.995 -40.120                       H
ATOM  12    C    PRO B  29      -2.203   0.437 -40.982                       C
ATOM  13    OC1  PRO B  29      -2.373   0.688 -42.152                       O
ATOM  14    HCO  PRO B  29      -2.227   1.208 -40.190                       H
ATOM  0     N    PRO C   1       1.990   3.037  44.208                       N
ATOM  1     H1   PRO C   1       1.167   2.413  44.242                       H
ATOM  2     CD   PRO C   1       3.159   2.370  44.868                       C
ATOM  3     HD1  PRO C   1       2.827   1.615  45.603                       H
ATOM  4     HD2  PRO C   1       3.710   3.155  45.429                       H
ATOM  5     CG   PRO C   1       4.045   1.763  43.763                       C
ATOM  6     HG1  PRO C   1       3.807   0.695  43.610                       H
ATOM  7     HG2  PRO C   1       5.114   1.814  44.019                       H
ATOM  8     CB   PRO C   1       3.721   2.571  42.497                       C
ATOM  9     HB1  PRO C   1       3.663   1.911  41.603                       H
ATOM  10    HB2  PRO C   1       4.509   3.311  42.284                       H
ATOM  11    CA   PRO C   1       2.371   3.283  42.779                       C
ATOM  12    HA   PRO C   1       2.465   4.390  42.665                       H
ATOM  13    C    PRO C   1       1.311   2.723  41.820                       C
ATOM  14    O    PRO C   1       0.476   1.899  42.174                       O
ATOM  0     N    PRO C   2       1.294   3.219  40.515                       N
ATOM  1     CD   PRO C   2       2.222   4.256  39.974                       C
ATOM  2     HD1  PRO C   2       3.278   3.920  40.045                       H
ATOM  3     HD2  PRO C   2       2.109   5.190  40.560                       H
ATOM  4     CG   PRO C   2       1.803   4.445  38.501                       C
ATOM  5     HG1  PRO C   2       2.670   4.307  37.831                       H
ATOM  6     HG2  PRO C   2       1.421   5.464  38.321                       H
ATOM  7     CB   PRO C   2       0.714   3.396  38.205                       C
ATOM  8     HB1  PRO C   2       0.896   2.872  37.241                       H
ATOM  9     HB2  PRO C   2      -0.284   3.865  38.081                       H
ATOM  10    CA   PRO C   2       0.705   2.430  39.403                       C
ATOM  11    HA   PRO C   2      -0.330   2.088  39.687                       H
ATOM  12    C    PRO C   2       1.594   1.192  39.120                       C
ATOM  13    O    PRO C   2       2.702   1.020  39.608                       O
ATOM  0     N    GLY C   3       1.023   0.275  38.249                       N
ATOM  1     H    GLY C   3       0.048   0.353  37.927                       H
ATOM  2     CA   GLY C   3       1.704  -0.980  37.916                       C
ATOM  3     HA1  GLY C   3       0.929  -1.791  37.872                       H
ATOM  4     HA2  GLY C   3       2.441  -1.247  38.719                       H
ATOM  5     C    GLY C   3       2.464  -0.884  36.592                       C
ATOM  6     O    GLY C   3       2.879   0.145  36.075                       O
ATOM  0     N    PRO C   4       2.773  -2.116  35.996                       N
ATOM  1     CD   PRO C   4       2.177  -3.426  36.389                       C
ATOM  2     HD1  PRO C   4       1.077  -3.328  36.529                       H
ATOM  3     HD2  PRO C   4       2.620  -3.748  37.356                       H
ATOM  4     CG   PRO C   4       2.520  -4.405  35.254                       C
ATOM  5     HG1  PRO C   4       1.625  -4.613  34.637                       H
ATOM  6     HG2  PRO C   4       2.865  -5.376  35.646                       H
ATOM  7     CB   PRO C   4       3.608  -3.720  34.409                       C
ATOM  8     HB1  PRO C   4       3.475  -3.953  33.335                       H
ATOM  9     HB2  PRO C   4       4.612  -4.079  34.693                       H
ATOM  10    CA   PRO C   4       3.482  -2.206  34.696                       C
ATOM  11    HA   PRO C   4       4.478  -1.703  34.806                       H
ATOM  12    C    PRO C   4       2.692  -1.532  33.553                       C
ATOM  13    O    PRO C   4       1.513  -1.199  33.636                       O
ATOM  0     N    PRO C   5       3.395  -1.288  32.381                       N
ATOM  1     CD   PRO C   5       4.827  -1.609  32.120                       C
ATOM  2     HD1  PRO C   5       4.971  -2.709  32.048                       H
ATOM  3     HD2  PRO C   5       5.460  -1.229  32.948                       H
ATOM  4     CG   PRO C   5       5.153  -0.918  30.783                       C
ATOM  5     HG1  PRO C   5       5.819  -1.542  30.163                       H
ATOM  6     HG2  PRO C   5       5.684   0.036  30.951                       H
ATOM  7     CB   PRO C   5       3.808  -0.672  30.077                       C
ATOM  8     HB1  PRO C   5       3.669  -1.386  29.242                       H
ATOM  9     HB2  PRO C   5       3.756   0.337  29.614                       H
ATOM  10    CA   PRO C   5       2.708  -0.838  31.144                       C
ATOM  11    HA   PRO C   5       2.186   0.144  31.367                       H
ATOM  12    C    PRO C   5       1.655  -1.873  30.687                       C
ATOM  13    O    PRO C   5       1.561  -3.016  31.114                       O
ATOM  0     N    GLY C   6       0.768  -1.356  29.756                       N
ATOM  1     H    GLY C   6       0.808  -0.376  29.454                       H
ATOM  2     CA   GLY C   6      -0.346  -2.148  29.231                       C
ATOM  3     HA1  GLY C   6      -1.221  -1.450  29.055                       H
ATOM  4     HA2  GLY C   6      -0.687  -2.894  29.994                       H
ATOM  5     C    GLY C   6      -0.001  -2.878  27.937                       C
ATOM  6     O    GLY C   6       1.118  -3.022  27.461                       O
ATOM  0     N    PRO C   7      -1.102  -3.461  27.300                       N
ATOM  1     CD   PRO C   7      -2.524  -3.252  27.706                       C
ATOM  2     HD1  PRO C   7      -2.717  -2.179  27.945                       H
ATOM  3     HD2  PRO C   7      -2.726  -3.832  28.633                       H
ATOM  4     CG   PRO C   7      -3.373  -3.746  26.525                       C
ATOM  5     HG1  PRO C   7      -3.780  -2.888  25.955                       H
ATOM  6     HG2  PRO C   7      -4.239  -4.335  26.869                       H
ATOM  7     CB   PRO C   7      -2.436  -4.585  25.639                       C
ATOM  8     HB1  PRO C   7      -2.671  -4.431  24.567                       H
ATOM  9     HB2  PRO C   7      -2.562  -5.663  25.843                       H
ATOM  10    CA   PRO C   7      -0.996  -4.140  25.987                       C
ATOM  11    HA   PRO C   7      -0.292  -5.006  26.094                       H
ATOM  12    C    PRO C   7      -0.481  -3.191  24.885                       C
ATOM  13    O    PRO C   7      -0.433  -1.968  24.998                       O
ATOM  0     N    PRO C   8      -0.036  -3.785  23.711                       N
ATOM  1     CD   PRO C   8      -0.001  -5.244  23.410                       C
ATOM  2     HD1  PRO C   8      -1.032  -5.649  23.323                       H
ATOM  3     HD2  PRO C   8       0.519  -5.787  24.225                       H
ATOM  4     CG   PRO C   8       0.753  -5.354  22.071                       C
ATOM  5     HG1  PRO C   8       0.314  -6.137  21.430                       H
ATOM  6     HG2  PRO C   8       1.808  -5.638  22.235                       H
ATOM  7     CB   PRO C   8       0.661  -3.972  21.401                       C
ATOM  8     HB1  PRO C   8      -0.054  -3.995  20.556                       H
ATOM  9     HB2  PRO C   8       1.630  -3.655  20.957                       H
ATOM  10    CA   PRO C   8       0.211  -2.978  22.488                       C
ATOM  11    HA   PRO C   8       1.013  -2.212  22.721                       H
ATOM  12    C    PRO C   8      -1.081  -2.246  22.053                       C
ATOM  13    O    PRO C   8      -2.204  -2.515  22.462                       O
ATOM  0     N    GLY C   9      -0.841  -1.209  21.169                       N
ATOM  1     H    GLY C   9       0.110  -0.947  20.883                       H
ATOM  2     CA   GLY C   9      -1.924  -0.370  20.652                       C
ATOM  3     HA1  GLY C   9      -1.512   0.673  20.485                       H
ATOM  4     HA2  GLY C   9      -2.736  -0.266  21.417                       H
ATOM  5     C    GLY C   9      -2.523  -0.905  19.356                       C
ATOM  6     O    GLY C   9      -2.332  -2.015  18.874                       O
ATOM  0     N    PRO C  10      -3.407  -0.021  18.726                       N
ATOM  1     CD   PRO C  10      -3.624   1.398  19.139                       C
ATOM  2     HD1  PRO C  10      -2.654   1.898  19.370                       H
ATOM  3     HD2  PRO C  10      -4.230   1.416  20.070                       H
ATOM  4     CG   PRO C  10      -4.357   2.066  17.966                       C
ATOM  5     HG1  PRO C  10      -3.664   2.713  17.394                       H
ATOM  6     HG2  PRO C  10      -5.175   2.716  18.319                       H
ATOM  7     CB   PRO C  10      -4.886   0.926  17.078                       C
ATOM  8     HB1  PRO C  10      -4.813   1.201  16.007                       H
ATOM  9     HB2  PRO C  10      -5.951   0.728  17.286                       H
ATOM  10    CA   PRO C  10      -4.034  -0.319  17.416                       C
ATOM  11    HA   PRO C  10      -4.652  -1.247  17.525                       H
ATOM  12    C    PRO C  10      -2.983  -0.528  16.306                       C
ATOM  13    O    PRO C  10      -1.799  -0.217  16.413                       O
ATOM  0     N    PRO C  11      -3.427  -1.122  15.132                       N
ATOM  1     CD   PRO C  11      -4.814  -1.581  14.837                       C
ATOM  2     HD1  PRO C  11      -5.505  -0.714  14.763                       H
ATOM  3     HD2  PRO C  11      -5.173  -2.246  15.649                       H
ATOM  4     CG   PRO C  11      -4.707  -2.322  13.490                       C
ATOM  5     HG1  PRO C  11      -5.590  -2.127  12.857                       H
ATOM  6     HG2  PRO C  11      -4.667  -3.415  13.644                       H
ATOM  7     CB   PRO C  11      -3.419  -1.821  12.814                       C
ATOM  8     HB1  PRO C  11      -3.658  -1.136  11.978                       H
ATOM  9     HB2  PRO C  11      -2.834  -2.649  12.358                       H
ATOM  10    CA   PRO C  11      -2.592  -1.108  13.902                       C
ATOM  11    HA   PRO C  11      -1.623  -1.651  14.123                       H
ATOM  12    C    PRO C  11      -2.278   0.346  13.479                       C
ATOM  13    O    PRO C  11      -2.863   1.335  13.901                       O
ATOM  0     N    GLY C  12      -1.221   0.432  12.589                       N
ATOM  1     H    GLY C  12      -0.693  -0.397  12.290                       H
ATOM  2     CA   GLY C  12      -0.746   1.719  12.079                       C
ATOM  3     HA1  GLY C  12       0.371   1.636  11.901                       H
ATOM  4     HA2  GLY C  12      -0.882   2.519  12.851                       H
ATOM  5     C    GLY C  12      -1.447   2.142  10.792                       C
ATOM  6     O    GLY C  12      -2.453   1.633  10.315                       O
ATOM  0     N    PRO C  13      -0.872   3.254  10.166                       N
ATOM  1     CD   PRO C  13       0.422   3.880  10.571                       C
ATOM  2     HD1  PRO C  13       1.191   3.101  10.786                       H
ATOM  3     HD2  PRO C  13       0.268   4.455  11.509                       H
ATOM  4     CG   PRO C  13       0.830   4.791   9.403                       C
ATOM  5     HG1  PRO C  13       1.648   4.328   8.818                       H
ATOM  6     HG2  PRO C  13       1.209   5.762   9.761                       H
ATOM  7     CB   PRO C  13      -0.425   4.963   8.528                       C
ATOM  8     HB1  PRO C  13      -0.151   4.985   7.454                       H
ATOM  9     HB2  PRO C  13      -0.931   5.918   8.749                       H
ATOM  10    CA   PRO C  13      -1.355   3.774   8.864                       C
ATOM  11    HA   PRO C  13      -2.425   4.086   8.985                       H
ATOM  12    C    PRO C  13      -1.251   2.719   7.744                       C
ATOM  13    O    PRO C  13      -0.601   1.680   7.837                       O
ATOM  0     N    PRO C  14      -1.960   2.976   6.578                       N
ATOM  1     CD   PRO C  14      -2.814   4.165   6.300                       C
ATOM  2     HD1  PRO C  14      -2.193   5.084   6.229                       H
ATOM  3     HD2  PRO C  14      -3.549   4.302   7.119                       H
ATOM  4     CG   PRO C  14      -3.499   3.854   4.956                       C
ATOM  5     HG1  PRO C  14      -3.582   4.760   4.332                       H
ATOM  6     HG2  PRO C  14      -4.529   3.488   5.114                       H
ATOM  7     CB   PRO C  14      -2.642   2.780   4.264                       C
ATOM  8     HB1  PRO C  14      -2.066   3.220   3.427                       H
ATOM  9     HB2  PRO C  14      -3.261   1.979   3.804                       H
ATOM  10    CA   PRO C  14      -1.707   2.194   5.339                       C
ATOM  11    HA   PRO C  14      -1.935   1.105   5.552                       H
ATOM  12    C    PRO C  14      -0.229   2.330   4.906                       C
ATOM  13    O    PRO C  14       0.544   3.180   5.330                       O
ATOM  0     N    GLY C  15       0.161   1.356   4.004                       N
ATOM  1     H    GLY C  15      -0.475   0.608   3.702                       H
ATOM  2     CA   GLY C  15       1.528   1.291   3.481                       C
ATOM  3     HA1  GLY C  15       1.779   0.203   3.291                       H
ATOM  4     HA2  GLY C  15       2.257   1.652   4.252                       H
ATOM  5     C    GLY C  15       1.712   2.099   2.201                       C
ATOM  6     O    GLY C  15       0.922   2.912   1.738                       O
ATOM  0     N    PRO C  16       2.940   1.889   1.563                       N
ATOM  1     CD   PRO C  16       3.927   0.838   1.951                       C
ATOM  2     HD1  PRO C  16       3.416  -0.131   2.159                       H
ATOM  3     HD2  PRO C  16       4.437   1.147   2.889                       H
ATOM  4     CG   PRO C  16       4.910   0.734   0.774                       C
ATOM  5     HG1  PRO C  16       4.709  -0.180   0.181                       H
ATOM  6     HG2  PRO C  16       5.952   0.659   1.124                       H
ATOM  7     CB   PRO C  16       4.693   1.991  -0.086                       C
ATOM  8     HB1  PRO C  16       4.789   1.748  -1.163                       H
ATOM  9     HB2  PRO C  16       5.454   2.758   0.139                       H
ATOM  10    CA   PRO C  16       3.282   2.519   0.266                       C
ATOM  11    HA   PRO C  16       3.260   3.632   0.398                       H
ATOM  12    C    PRO C  16       2.298   2.114  -0.852                       C
ATOM  13    O    PRO C  16       1.502   1.182  -0.761                       O
ATOM  0     N    PRO C  17       2.323   2.878  -2.011                       N
ATOM  1     CD   PRO C  17       3.200   4.052  -2.284                       C
ATOM  2     HD1  PRO C  17       4.262   3.735  -2.367                       H
ATOM  3     HD2  PRO C  17       3.118   4.786  -1.457                       H
ATOM  4     CG   PRO C  17       2.687   4.627  -3.618                       C
ATOM  5     HG1  PRO C  17       3.522   4.984  -4.245                       H
ATOM  6     HG2  PRO C  17       2.030   5.498  -3.446                       H
ATOM  7     CB   PRO C  17       1.914   3.494  -4.316                       C
ATOM  8     HB1  PRO C  17       2.500   3.084  -5.161                       H
ATOM  9     HB2  PRO C  17       0.961   3.850  -4.765                       H
ATOM  10    CA   PRO C  17       1.642   2.417  -3.249                       C
ATOM  11    HA   PRO C  17       0.535   2.308  -3.028                       H
ATOM  12    C    PRO C  17       2.208   1.050  -3.700                       C
ATOM  13    O    PRO C  17       3.252   0.560  -3.289                       O
ATOM  0     N    GLY C  18       1.386   0.398  -4.602                       N
ATOM  1     H    GLY C  18       0.480   0.785  -4.892                       H
ATOM  2     CA   GLY C  18       1.724  -0.922  -5.139                       C
ATOM  3     HA1  GLY C  18       0.757  -1.483  -5.328                       H
ATOM  4     HA2  GLY C  18       2.288  -1.520  -4.377                       H
ATOM  5     C    GLY C  18       2.543  -0.849  -6.423                       C
ATOM  6     O    GLY C  18       3.081   0.150  -6.883                       O
ATOM  0     N    PRO C  19       2.703  -2.080  -7.071                       N
ATOM  1     CD   PRO C  19       1.996  -3.338  -6.687                       C
ATOM  2     HD1  PRO C  19       0.920  -3.140  -6.470                       H
ATOM  3     HD2  PRO C  19       2.449  -3.738  -5.755                       H
ATOM  4     CG   PRO C  19       2.181  -4.299  -7.871                       C
ATOM  5     HG1  PRO C  19       1.246  -4.374  -8.460                       H
ATOM  6     HG2  PRO C  19       2.421  -5.319  -7.530                       H
ATOM  7     CB   PRO C  19       3.313  -3.712  -8.734                       C
ATOM  8     HB1  PRO C  19       3.105  -3.870  -9.811                       H
ATOM  9     HB2  PRO C  19       4.273  -4.211  -8.516                       H
ATOM  10    CA   PRO C  19       3.397  -2.210  -8.374                       C
ATOM  11    HA   PRO C  19       4.454  -1.857  -8.247                       H
ATOM  12    C    PRO C  19       2.710  -1.385  -9.482                       C
ATOM  13    O    PRO C  19       1.584  -0.903  -9.380                       O
ATOM  0     N    PRO C  20       3.438  -1.173 -10.645                       N
ATOM  1     CD   PRO C  20       4.818  -1.658 -10.931                       C
ATOM  2     HD1  PRO C  20       4.831  -2.765 -11.025                       H
ATOM  3     HD2  PRO C  20       5.500  -1.370 -10.105                       H
ATOM  4     CG   PRO C  20       5.207  -0.984 -12.260                       C
ATOM  5     HG1  PRO C  20       5.793  -1.667 -12.898                       H
ATOM  6     HG2  PRO C  20       5.843  -0.099 -12.082                       H
ATOM  7     CB   PRO C  20       3.891  -0.576 -12.946                       C
ATOM  8     HB1  PRO C  20       3.668  -1.252 -13.795                       H
ATOM  9     HB2  PRO C  20       3.944   0.442 -13.388                       H
ATOM  10    CA   PRO C  20       2.788  -0.645 -11.872                       C
ATOM  11    HA   PRO C  20       2.360   0.377 -11.638                       H
ATOM  12    C    PRO C  20       1.646  -1.583 -12.326                       C
ATOM  13    O    PRO C  20       1.486  -2.730 -11.926                       O
ATOM  0     N    GLY C  21       0.775  -0.980 -13.217                       N
ATOM  1     H    GLY C  21       0.878   0.004 -13.495                       H
ATOM  2     CA   GLY C  21      -0.396  -1.682 -13.747                       C
ATOM  3     HA1  GLY C  21      -1.213  -0.917 -13.925                       H
ATOM  4     HA2  GLY C  21      -0.798  -2.399 -12.986                       H
ATOM  5     C    GLY C  21      -0.102  -2.436 -15.040                       C
ATOM  6     O    GLY C  21       1.005  -2.655 -15.516                       O
ATOM  0     N    PRO C  22      -1.240  -2.944 -15.677                       N
ATOM  1     CD   PRO C  22      -2.645  -2.640 -15.272                       C
ATOM  2     HD1  PRO C  22      -2.769  -1.555 -15.046                       H
ATOM  3     HD2  PRO C  22      -2.883  -3.197 -14.340                       H
ATOM  4     CG   PRO C  22      -3.525  -3.092 -16.448                       C
ATOM  5     HG1  PRO C  22      -3.876  -2.216 -17.028                       H
ATOM  6     HG2  PRO C  22      -4.426  -3.620 -16.099                       H
ATOM  7     CB   PRO C  22      -2.643  -3.998 -17.326                       C
ATOM  8     HB1  PRO C  22      -2.870  -3.844 -18.399                       H
ATOM  9     HB2  PRO C  22      -2.835  -5.063 -17.107                       H
ATOM  10    CA   PRO C  22      -1.178  -3.639 -16.984                       C
ATOM  11    HA   PRO C  22      -0.528  -4.545 -16.871                       H
ATOM  12    C    PRO C  22      -0.605  -2.732 -18.094                       C
ATOM  13    O    PRO C  22      -0.474  -1.515 -17.987                       O
ATOM  0     N    PRO C  23      -0.205  -3.360 -19.266                       N
ATOM  1     CD   PRO C  23      -0.265  -4.819 -19.560                       C
ATOM  2     HD1  PRO C  23      -1.321  -5.155 -19.652                       H
ATOM  3     HD2  PRO C  23       0.213  -5.392 -18.740                       H
ATOM  4     CG   PRO C  23       0.487  -4.985 -20.894                       C
ATOM  5     HG1  PRO C  23       0.004  -5.743 -21.533                       H
ATOM  6     HG2  PRO C  23       1.521  -5.334 -20.722                       H
ATOM  7     CB   PRO C  23       0.487  -3.604 -21.572                       C
ATOM  8     HB1  PRO C  23      -0.231  -3.583 -22.414                       H
ATOM  9     HB2  PRO C  23       1.473  -3.355 -22.021                       H
ATOM  10    CA   PRO C  23       0.109  -2.576 -20.489                       C
ATOM  11    HA   PRO C  23       0.968  -1.875 -20.252                       H
ATOM  12    C    PRO C  23      -1.118  -1.745 -20.931                       C
ATOM  13    O    PRO C  23      -2.261  -1.920 -20.530                       O
ATOM  0     N    GLY C  24      -0.789  -0.732 -21.817                       N
ATOM  1     H    GLY C  24       0.183  -0.546 -22.092                       H
ATOM  2     CA   GLY C  24      -1.794   0.205 -22.323                       C
ATOM  3     HA1  GLY C  24      -1.289   1.205 -22.494                       H
ATOM  4     HA2  GLY C  24      -2.586   0.383 -21.551                       H
ATOM  5     C    GLY C  24      -2.451  -0.267 -23.616                       C
ATOM  6     O    GLY C  24      -2.357  -1.381 -24.116                       O
ATOM  0     N    PRO C  25      -3.265   0.695 -24.226                       N
ATOM  1     CD   PRO C  25      -3.369   2.120 -23.791                       C
ATOM  2     HD1  PRO C  25      -2.362   2.543 -23.566                       H
ATOM  3     HD2  PRO C  25      -3.960   2.171 -22.851                       H
ATOM  4     CG   PRO C  25      -4.066   2.858 -24.944                       C
ATOM  5     HG1  PRO C  25      -3.333   3.455 -25.522                       H
ATOM  6     HG2  PRO C  25      -4.824   3.566 -24.571                       H
ATOM  7     CB   PRO C  25      -4.697   1.772 -25.834                       C
ATOM  8     HB1  PRO C  25      -4.633   2.054 -26.903                       H
ATOM  9     HB2  PRO C  25      -5.768   1.647 -25.599                       H
ATOM  10    CA   PRO C  25      -3.929   0.464 -25.530                       C
ATOM  11    HA   PRO C  25      -4.610  -0.421 -25.428                       H
ATOM  12    C    PRO C  25      -2.907   0.195 -26.656                       C
ATOM  13    O    PRO C  25      -1.701   0.411 -26.556                       O
ATOM  0     N    PRO C  26      -3.406  -0.338 -27.836                       N
ATOM  1     CD   PRO C  26      -4.824  -0.692 -28.124                       C
ATOM  2     HD1  PRO C  26      -5.439   0.228 -28.229                       H
ATOM  3     HD2  PRO C  26      -5.243  -1.297 -27.295                       H
ATOM  4     CG   PRO C  26      -4.774  -1.481 -29.446                       C
ATOM  5     HG1  PRO C  26      -5.633  -1.233 -30.092                       H
ATOM  6     HG2  PRO C  26      -4.830  -2.569 -29.256                       H
ATOM  7     CB   PRO C  26      -3.442  -1.117 -30.125                       C
ATOM  8     HB1  PRO C  26      -3.615  -0.444 -30.986                       H
ATOM  9     HB2  PRO C  26      -2.929  -2.011 -30.542                       H
ATOM  10    CA   PRO C  26      -2.564  -0.440 -29.055                       C
ATOM  11    HA   PRO C  26      -1.649  -1.062 -28.805                       H
ATOM  12    C    PRO C  26      -2.109   0.961 -29.525                       C
ATOM  13    O    PRO C  26      -2.585   2.019 -29.134                       O
ATOM  0     N    GLY C  27      -1.063   0.913 -30.432                       N
ATOM  1     H    GLY C  27      -0.614   0.028 -30.699                       H
ATOM  2     CA   GLY C  27      -0.443   2.131 -30.956                       C
ATOM  3     HA1  GLY C  27       0.645   1.909 -31.183                       H
ATOM  4     HA2  GLY C  27      -0.433   2.929 -30.168                       H
ATOM  5     C    GLY C  27      -1.131   2.662 -32.206                       C
ATOM  6     O    GLY C  27      -2.144   2.218 -32.736                       O
ATOM  0     N    PRO C  28      -0.525   3.794 -32.763                       N
ATOM  1     CD   PRO C  28       0.778   4.385 -32.338                       C
ATOM  2     HD1  PRO C  28       1.577   3.610 -32.336                       H
ATOM  3     HD2  PRO C  28       0.693   4.770 -31.298                       H
ATOM  4     CG   PRO C  28       1.065   5.512 -33.345                       C
ATOM  5     HG1  PRO C  28       1.845   5.199 -34.064                       H
ATOM  6     HG2  PRO C  28       1.446   6.416 -32.840                       H
ATOM  7     CB   PRO C  28      -0.263   5.792 -34.075                       C
ATOM  8     HB1  PRO C  28      -0.093   6.058 -35.132                       H
ATOM  9     HB2  PRO C  28      -0.787   6.649 -33.614                       H
ATOM  10    CA   PRO C  28      -1.109   4.509 -33.921                       C
ATOM  11    HA   PRO C  28      -2.183   4.734 -33.686                       H
ATOM  12    C    PRO C  28      -1.017   3.655 -35.204                       C
ATOM  13    O    PRO C  28      -0.097   2.877 -35.442                       O
ATOM  0     N    PRO C  29      -2.004   3.841 -36.162                       N
ATOM  1     CD   PRO C  29      -3.162   4.775 -36.051                       C
ATOM  2     HD1  PRO C  29      -2.797   5.817 -35.939                       H
ATOM  3     HD2  PRO C  29      -3.772   4.525 -35.155                       H
ATOM  4     CG   PRO C  29      -3.970   4.592 -37.347                       C
ATOM  5     HG1  PRO C  29      -3.786   5.430 -38.045                       H
ATOM  6     HG2  PRO C  29      -5.055   4.586 -37.150                       H
ATOM  7     CB   PRO C  29      -3.500   3.265 -37.965                       C
ATOM  8     HB1  PRO C  29      -3.435   3.328 -39.070                       H
ATOM  9     HB2  PRO C  29      -4.217   2.445 -37.738                       H
ATOM  10    CA   PRO C  29      -2.125   2.951 -37.352                       C
ATOM  11    HA   PRO C  29      -2.060   1.862 -37.030                       H
ATOM  12    C    PRO C  29      -0.964   3.175 -38.327                       C
ATOM  13    OC1  PRO C  29      -1.114   3.048 -39.523                       O
ATOM  14    HCO  PRO C  29       0.009   3.416 -37.871                       H
"""

@attr(speed = 'slow' )
class OptimizedPpgTestCase( unittest.TestCase ):

    def setUp(self):
        self.sys = System.from_pdb_string( _ppg_string )
        self.ch1 = self.sys[0]
        self.ch1.connect_residues()

        ch1 = self.ch1
        rel_resid = [1, 2]
        rel_resname = ['PRO']
        relevant = [res for res in ch1.molecules if res.res_id in rel_resid and res.res_name in rel_resname]

        for res in relevant:
            res.gather_ready( level = 2, residue = 1 )
            res.gather_ready( level = 2, concap = 1 )

        res1 = ch1[0].ready
        res2 = ch1[1].ready
        con1 = ch1[0].concap
        for each in [res1, res2, con1]:
            each.populate_bonds()
        self.res1  = res1 
        self.res2  = res2 
        self.con1 = con1
        self.res = ch1[0]

        self.res.to_AU()
        self.res1.to_AU()
        self.res2.to_AU()
        self.con1.to_AU()


    def test_atoms_and_bonds(self):
        res1  = self.res1
        res2  = self.res2  
        con1  = self.con1 

#bonds is actually counted twice, the 'Atom.get_ats_and_bonds' method returns
#unique bonds only not duplicates
        assert len(res1) == 21
        assert len(res1.bonds) == 42

        assert len(res2) == 26
        assert len(res2.bonds) == 52

        assert len(con1) == 12
        assert len(con1.bonds) == 22

    def test_water_attach_from_targz(self):
        w1 = Water.get_standard()
        w2 = Water.get_standard()

        w1.populate_bonds()
        w1.props_from_targz( WATER_FILE_TARGZ, bonds = True, maxl = 1, pol = 22, hyp = 1 )
        w2.props_from_targz( WATER_FILE_TARGZ, bonds = False, maxl = 1, pol = 22, hyp = 1 )

        np.testing.assert_allclose( w1.p.q, w2.p.q, atol = 1e-7 )
        np.testing.assert_allclose( w1.p.d, w2.p.d, atol = 1e-7 )
        np.testing.assert_allclose( w1.p.a, w2.p.a, atol = 1e-7 )
        np.testing.assert_allclose( w1.p.b, w2.p.b, atol = 1e-7 )

    def test_first_3_attach_from_targz(self):

        basis = ['6_31pgp_loprop']
        sys = System.from_pdb_string( _ppg_string )
        sys.connect_residues()

        rel_resid = [1, 2]
        rel_resname = ['PRO']

        relevant = [res for res in sys.molecules if res.res_id in rel_resid and res.res_name in rel_resname]

        for res in relevant:
            res.gather_ready( level = 2, residue = 1 )
            res.gather_ready( level = 2, concap = 1 )

        res = sys[0][0]
        res1 = sys[0][0].ready
        res2 = sys[0][1].ready
        con1 = sys[0][0].concap

        res1.props_from_targz( PRO_1_FILE_TARGZ, maxl = 1, bonds = 1 )
        res2.props_from_targz( PRO_2_FILE_TARGZ, maxl = 1, bonds = 1 )
        con1.props_from_targz( CON_1_FILE_TARGZ, maxl = 1, bonds = 1 )

        n = res1.get_atom_by_pdbname( 'N', dup = 1 )[0]
        cg = res1.CG
        res.mfcc_props()
        np.testing.assert_allclose( res.N.p.q, n.q )
        np.testing.assert_allclose( res.CG.p.q, cg.q )
#hand calc. a_zz bond midpoint of C-N all atom contr. C first proline
        np.testing.assert_allclose( res.C.p.a[5], 8.0552875, atol=1e-7 )


if __name__ == '__main__':
    unittest.main()