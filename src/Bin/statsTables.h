/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    statsTables.h
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Pre-computed critical values for t-test  
      Check:   README and gpl.txt
      Copyright (C) 2008 Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
*/

/*************************************************************************************

    Program's name: els-ptsp

    Estimation-based Local Search for the PTSP 

	Copyright (C) 2008 Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Email: {pbalapra,mbiro,stuetzle}@ulb.ac.be
    Mail address:	Prasanna Balaprakash, 
					IRIDIA, Universite Libre de Bruxelles 
					50, Av. F. Roosevelt, CP 194/6 
					B-1050 Brussels, Belgium 
					http://iridia.ulb.ac.be/~prasanna

*****************************************************************************************/

double
percentage_points_t_distribution[122][4]=
  {
    {-1.000, -1.000, -1.000, -1.000},
    {6.314, 12.706, 31.821, 63.657},
    {2.920, 4.303, 6.965, 9.925},
    {2.353, 3.182, 4.541, 5.841},
    {2.132, 2.776, 3.747, 4.604},
    {2.015, 2.571, 3.365, 4.032},
    {1.943, 2.447, 3.143, 3.707},
    {1.895, 2.365, 2.998, 3.499},
    {1.860, 2.306, 2.896, 3.355},
    {1.833, 2.262, 2.821, 3.250},
    {1.812, 2.228, 2.764, 3.169},
    {1.796, 2.201, 2.718, 3.106},
    {1.782, 2.179, 2.681, 3.055},
    {1.771, 2.160, 2.650, 3.012},
    {1.761, 2.145, 2.624, 2.977},
    {1.753, 2.131, 2.602, 2.947},
    {1.746, 2.120, 2.583, 2.921},
    {1.740, 2.110, 2.567, 2.898},
    {1.734, 2.101, 2.552, 2.878},
    {1.729, 2.093, 2.539, 2.861},
    {1.725, 2.086, 2.528, 2.845},
    {1.721, 2.080, 2.518, 2.831},
    {1.717, 2.074, 2.508, 2.819},
    {1.714, 2.069, 2.500, 2.807},
    {1.711, 2.064, 2.492, 2.797},
    {1.708, 2.060, 2.485, 2.787},
    {1.706, 2.056, 2.479, 2.779},
    {1.703, 2.052, 2.473, 2.771},
    {1.701, 2.048, 2.467, 2.763},
    {1.699, 2.045, 2.462, 2.756},
    {1.697, 2.042, 2.457, 2.750},
    {1.696, 2.040, 2.453, 2.744},
    {1.694, 2.037, 2.449, 2.738},
    {1.692, 2.035, 2.445, 2.733},
    {1.691, 2.032, 2.441, 2.728},
    {1.690, 2.030, 2.438, 2.724},
    {1.688, 2.028, 2.434, 2.719},
    {1.687, 2.026, 2.431, 2.715},
    {1.686, 2.024, 2.429, 2.712},
    {1.685, 2.023, 2.426, 2.708},
    {1.684, 2.021, 2.423, 2.704},
    {1.683, 2.020, 2.421, 2.701},
    {1.682, 2.018, 2.418, 2.698},
    {1.681, 2.017, 2.416, 2.695},
    {1.680, 2.015, 2.414, 2.692},
    {1.679, 2.014, 2.412, 2.690},
    {1.679, 2.013, 2.410, 2.687},
    {1.678, 2.012, 2.408, 2.685},
    {1.677, 2.011, 2.407, 2.682},
    {1.677, 2.010, 2.405, 2.680},
    {1.676, 2.009, 2.403, 2.678},
    {1.675, 2.008, 2.402, 2.676},
    {1.675, 2.007, 2.400, 2.674},
    {1.674, 2.006, 2.399, 2.672},
    {1.674, 2.005, 2.397, 2.670},
    {1.673, 2.004, 2.396, 2.668},
    {1.673, 2.003, 2.395, 2.667},
    {1.672, 2.002, 2.394, 2.665},
    {1.672, 2.002, 2.392, 2.663},
    {1.671, 2.001, 2.391, 2.662},
    {1.671, 2.000, 2.390, 2.660},
    {1.670, 2.000, 2.389, 2.659},
    {1.670, 1.999, 2.388, 2.657},
    {1.669, 1.998, 2.387, 2.656},
    {1.669, 1.998, 2.386, 2.655},
    {1.669, 1.997, 2.385, 2.654},
    {1.668, 1.997, 2.384, 2.652},
    {1.668, 1.996, 2.383, 2.651},
    {1.668, 1.995, 2.382, 2.650},
    {1.667, 1.995, 2.382, 2.649},
    {1.667, 1.994, 2.381, 2.648},
    {1.667, 1.994, 2.380, 2.647},
    {1.666, 1.993, 2.379, 2.646},
    {1.666, 1.993, 2.379, 2.645},
    {1.666, 1.993, 2.378, 2.644},
    {1.665, 1.992, 2.377, 2.643},
    {1.665, 1.992, 2.376, 2.642},
    {1.665, 1.991, 2.376, 2.641},
    {1.665, 1.991, 2.375, 2.640},
    {1.664, 1.990, 2.374, 2.640},
    {1.664, 1.990, 2.374, 2.639},
    {1.664, 1.990, 2.373, 2.638},
    {1.664, 1.989, 2.373, 2.637},
    {1.663, 1.989, 2.372, 2.636},
    {1.663, 1.989, 2.372, 2.636},
    {1.663, 1.988, 2.371, 2.635},
    {1.663, 1.988, 2.370, 2.634},
    {1.663, 1.988, 2.370, 2.634},
    {1.662, 1.987, 2.369, 2.633},
    {1.662, 1.987, 2.369, 2.632},
    {1.662, 1.987, 2.368, 2.632},
    {1.662, 1.986, 2.368, 2.631},
    {1.662, 1.986, 2.368, 2.630},
    {1.661, 1.986, 2.367, 2.630},
    {1.661, 1.986, 2.367, 2.629},
    {1.661, 1.985, 2.366, 2.629},
    {1.661, 1.985, 2.366, 2.628},
    {1.661, 1.985, 2.365, 2.627},
    {1.661, 1.984, 2.365, 2.627},
    {1.660, 1.984, 2.365, 2.626},
    {1.660, 1.984, 2.364, 2.626},
    {1.660, 1.984, 2.364, 2.625},
    {1.660, 1.983, 2.363, 2.625},
    {1.660, 1.983, 2.363, 2.624},
    {1.660, 1.983, 2.363, 2.624},
    {1.659, 1.983, 2.362, 2.623},
    {1.659, 1.983, 2.362, 2.623},
    {1.659, 1.982, 2.362, 2.623},
    {1.659, 1.982, 2.361, 2.622},
    {1.659, 1.982, 2.361, 2.622},
    {1.659, 1.982, 2.361, 2.621},
    {1.659, 1.982, 2.360, 2.621},
    {1.659, 1.981, 2.360, 2.620},
    {1.658, 1.981, 2.360, 2.620},
    {1.658, 1.981, 2.360, 2.620},
    {1.658, 1.981, 2.359, 2.619},
    {1.658, 1.981, 2.359, 2.619},
    {1.658, 1.980, 2.359, 2.619},
    {1.658, 1.980, 2.358, 2.618},
    {1.658, 1.980, 2.358, 2.618},
    {1.658, 1.980, 2.358, 2.617},
    {1.645, 1.960, 2.326, 2.576}
  };

