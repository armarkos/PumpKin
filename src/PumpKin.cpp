//
//  PumpKin.cpp
//  PumpKin
//
//  Version 1.4
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


#include <iostream>
#include "pk_IO.h"
#include "pk_DataTypes.h"
#include "pk_Pathways.h"
#include "pk_Print.h"

using namespace std;

size_t n_R;         // number of reactions
size_t n_S;         // number of species
size_t n_t;         // number of time-steps
double t_init;      // initial time
double t_end;       // final time
double f_min;
double tau_lifetime;
int interest;
int max_path;
int max_bp;
bool global_kin;

int main(int argc, const char *argv[])
{
    // Reads the location of input folder
    string folder = "";
    if (argc > 1) folder = string(argv[1]);

    cout.precision(12);

    Print_license();

    In_data kinetics;               // Keeps all input data
    Rates   rates;                  // Keeps all total rates

    Read_kin(kinetics, folder);     // Reads input file
    Average_kin(kinetics, rates);  // Average input data

    vector<pathway> PATHS(n_R);     // The set of all pathways
    Initialize_PATHS(PATHS, kinetics, rates);
    Initialize_rates(PATHS, kinetics, rates);

    // All brenching points so far
    vector<int> all_brenching_points;

    // Run baby, run!
    Run(PATHS, all_brenching_points, kinetics, rates);

    Print_Results(PATHS, rates, kinetics, folder, all_brenching_points);

    return 0;
}














