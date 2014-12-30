//
//  pk_IO.h
//  PumpKin
//
//  Version 1.2
//
//  Created by Aram H. Markosyan on 9/21/13.
//  Copyright (c) 2013 - 2015 Aram H. Markosyan. All rights reserved.
//
//
// This file is part of PumpKin (see http://www.pumpkin-tool.org).
// Please see the Copyright_and_License file for the copyright notice,
// disclaimer, contact information and the License.
//
// PumpKin is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, contact Aram H. Markosyan at
// aram.math@gmail.com or write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

#ifndef __PumpKin__pk_IO__
#define __PumpKin__pk_IO__

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <dirent.h>
#include <algorithm>
#include "stdlib.h"
#include <iomanip>
#include "pk_Array.h"
#include "pk_DataTypes.h"


// Reads input kinetic model
void Read_kin(In_data      &kinetics,
              const string &folder);

// Reads the species from files
void Read_species(In_data      &kinetics,
                  const string &m_folder);

// Reads the stoichiometric matrix from files
void Read_matrix(In_data      &kinetics,
                 const string &m_folder);

// Checking if all values are numbers
bool Check_file_line(string &line);

// Reads the time steps from files
void Read_time(In_data      &kinetics,
               const string &m_folder);

// Reads the density from files
void Read_density(In_data      &kinetics,
                  const string &m_folder);

// Reads the rates from files
void Read_rates(In_data      &kinetics,
                const string &m_folder);

// Reads the the reactions from files
void Read_reactions(In_data      &kinetics,
                    const string &m_folder);

// Average input kinetic model
void Average_kin(In_data &kinetics,
                 Rates   &rates);

// Interpolation function
void interp_1(doublearray1d               &interpol,
              const double                &t,
              const doublearray1d         &T,
              const vector<doublearray1d> &Y);

// Print license information
void Print_license();



#endif /* defined(__PumpKin__pk_IO__) */





