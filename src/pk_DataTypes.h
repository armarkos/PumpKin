//
//  pk_DataTypes.h
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

#ifndef PumpKin_pk_DataTypes_h
#define PumpKin_pk_DataTypes_h

#include <string>
#include <vector>
#include "pk_Array.h"

using namespace std;

/*
 Expected files in input folder:

 1. qt_reactions_list.txt
 2. qt_species_list.txt
 3. qt_matrix.txt
 4. qt_conditions.txt
 5. qt_densities.txt
 6. qt_rates.txt
 7. input.txt
 */


extern size_t n_R;          // number of reactions
extern size_t n_S;          // number of species
extern size_t n_t;          // number of time-steps
extern double t_init;
extern double t_end;
extern double f_min;
extern double tau_lifetime;
extern int interest;
extern int max_path;
extern int max_bp;
extern bool global_kin;

// Keeps all input data
struct In_data
{
    vector<string>         reactions;
    vector<string>         species;
    doublearray2d          matrix;
    doublearray1d          time;
    vector<doublearray1d>  densities;   // vector<time>
    vector<doublearray1d>  rates;       // vector<time>
};

// Keeps total rates existing and deleted ones
struct Rates
{
    doublearray1d r_j;          // average rates
    doublearray1d c_i;          // average densities
    doublearray1d delta_c;      // change of conctentrations

    doublearray1d p_i;          // total rate of production
    doublearray1d d_i;          // total rate of destruction

    doublearray1d delta_i;      // the mean rate of concentration change
    doublearray1d D_i;          // auxiliary variable

    doublearray1d tilda_r_j;    // part of the rate of reactions deleted pathways
    doublearray1d tilda_p_i;    // rate of production of species deleted pathways
    doublearray1d tilda_d_i;    // rate of destruction of species deleted pathways

    vector<string> reactions;
    vector<string> species;
};


// Keeps pathways
struct pathway
{
    // Indices of reaction inside of path or "reaction sequence"
    vector<int> path;

    // Multiplicity of reactions from path
    vector<int> x_j;

    // Number of molecules of S_i producted/consumed
    vector<int> m_i;

    // Number of 'virtual' molecules produced/consumed.
    // This is used to implement the addition of pseudo-pathways in the
    // splitting into elementary pathways.
    vector<int> virtualm_i;

    // pathway rate
    double f_k;

    // deleted
    bool deleted;

    // is elementary
    bool simple;
};

// Keeps all elementary pathways
struct elementary
{
    vector<pathway> EL_pathways;
};

#endif
