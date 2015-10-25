//
//  pk_Pathways.h
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

#ifndef __PumpKin__pk_Pathways__
#define __PumpKin__pk_Pathways__

#include <iostream>
#include "pk_DataTypes.h"
#include "pk_Print.h"
#include <cmath>
#include <numeric>
#include <algorithm>
#include "glpk.h"


// Initialization of pathways
void Initialize_PATHS(vector<pathway> &PATHS,
                      const In_data   &kinetics,
                      const Rates     &rates);


// Calculate rates
void Initialize_rates(vector<pathway> &PATHS,
                      const In_data   &kinetics,
                      Rates           &rates);


// Find and sort brenching points
void Branching_points(vector<int> &bp,
                      const Rates &rates);

// Get the next brenching point, returns false if there is no available one
bool Get_bp(vector<int> &all_bp,
            const Rates &rates,
            int         &bp);

// Run baby, run!
void Run(vector<pathway> &PATHS,
         vector<int>     &all_brenching_points,
         const In_data   &kinetics,
         Rates           &rates);


// Merge pathways at this brenching point
void Merge_pathways(vector<pathway> &PATHS,
                    const int       &brenching_point,
                    const In_data   &kinetics,
                    Rates           &rates);


// Merge 2 pathways
void Merge_2_pathways(vector<pathway> &PATHS,
                      const pathway   &Pathway_producing,
                      const pathway   &Pathway_consuming,
                      const int       &brenching_point,
                      const In_data   &kinetics,
                      Rates           &rates);

// Update the p, d
void Calculate_rates(vector<pathway> &PATHS,
                     const In_data   &kinetics,
                     Rates           &rates);


// Sort and delete repeted reactions in pathway paths
void Sort_and_delete_pathway(vector<int> &path,
                             vector<int> &xj);


// Elimination of insignificant pathways
void Insignificant(vector<pathway> &PATHS,
                   Rates           &rates,
                   const In_data   &input);

// Updates the deleted rates after merging two pathways
void Update_deleted_after_merging(vector<pathway> &PATHS,
                                  pathway         &producing,
                                  pathway         &consuming,
                                  Rates           &rates,
                                  double           dtilde,
                                  double           ptilde,
                                  int              branching_point);

// Updates the deleted rates after going through a branching point
void Update_deleted_for_branching_point(vector<pathway> &PATHS,
                                        Rates           &rates,
                                        int              branching_point);

// Updates the contributions from deleted pathways, applying (43)-(45).
void Update_deleted(pathway         &p,
                    Rates           &rates,
                    double           f);

// Create elementary pathways
void Make_elementary(vector<pathway>    &PATHS,
                     vector<elementary> &Elementary,
                     const vector<int>  &all_brenching_points,
                     const Rates        &rates,
                     const In_data      &kinetics);


// Splitting into Sub-Pathways
void Splitting(vector<pathway>    &PATHS,
               vector<elementary> &Elementary,
               const vector<int>  &all_brenching_points,
               const Rates        &rates,
               const In_data      &kinetics);


// Implements step_1 from 5.5.1
void Step_1(const pathway     &P_n,
            elementary        &elementary_paths,
            const Rates       &rates,
            const In_data     &kinetics);

// Adds pseudo-pathways to elementaries for all branching points so far
void Add_pseudo_pathways(const pathway     &P_n,
                         elementary        &elementaries,
                         const vector<int> &all_brenching_points);

void Remove_pseudo_pathways(elementary &elementaries);

// Implements step_2 from 5.5.1
void Step_2(elementary        &e_path,
            const vector<int> &all_brenching_points,
            const In_data     &kinetics);


// Implements step_2_3 from 5.5.1
void Step_2_3(const elementary &e_path,
              elementary       &P_tilda,
              int               bp,
              const In_data    &kinetics);

// Combine 2 pathways
void Combine_2_pathways(pathway           &combination,
                        const pathway     &Pathway_producing,
                        const pathway     &Pathway_consuming,
                        const int         &brenching_point,
                        const In_data     &input_raw,
                        const elementary  &e_path);

// Checks the validity of formula in 2.3
void Condition_2_3(const elementary &e_path,
                   elementary       &P_tilda,
                   const pathway    &combine,
                   const pathway     &Pathway_producing,
                   const pathway     &Pathway_consuming,
                   int k, int l,
                   const In_data    &kinetics);


// Calculating w_k
void Get_w_k(pathway       &path,
             elementary    &element,
             doublearray1d &w_k);


// Creates square LP problem
void Creat_LP(const pathway    &pw,
              const elementary &element,
              doublearray2d    &A,
              doublearray1d    &b);


/*
 solve LP problem, where A is square matrix:

 z = \sum_i c_i * x_i -> min
 Ax = b
 x >= 0

 */
void Solve_LP(double        &z,
              doublearray1d &x,
              doublearray2d &A,
              doublearray1d &b,
              doublearray1d &c);

// GCD of x_j with f_k update
void Factorise (pathway &pw);

// Sort by rate
void Sort_PATH_by_rate(vector<pathway> &PATHS);

bool compare_pathway_rate(const pathway &a,
                          const pathway &b);

// Copy pathway
void Copy_pathway(const pathway &from,
                  pathway       &to);

// delete pathway
void Del_pathway(pathway &del);

// Checks if one \in in_two + in_three
bool Contains(const vector<int> &one,
              const vector<int> &in_two,
              const vector<int> &in_three);

// Checking if pathway_1 = pathway_2
bool Paths_are_equal(pathway pathway_1,
                     pathway pathway_2);

// Returns GCD
int GCD (vector<int> x_j);

int gcd_local(int a,
              int b);



#endif /* defined(__PumpKin__pk_Pathways__) */





