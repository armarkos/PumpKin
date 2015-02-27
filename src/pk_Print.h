//
//  pk_Print.h
//  PumpKin
//
//  Version 1.3
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

#ifndef __PumpKin__Print__
#define __PumpKin__Print__

#include <iostream>
#include "pk_DataTypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

// Printing individual pathway
void Print_pathway(const pathway &PATHWAY,
                   const In_data &input);


// Printing all pathways
void Print_PATHS(const vector<pathway> &PATHS,
                 const Rates           &rates);


// Printing double array in 2D
void Print_doublearray2D(const doublearray2d &A);


// Printing final report
void Print_report(const vector<pathway> &PATHS,
                  const vector<int>     &all_brenching_points);

// Printing all pathways according to specie of interest
void Print_PATHS_for_specie(const vector<pathway> &PATHS,
                            const Rates           &rates);

// Printing all reduced system
void Print_PATHS_reduce(const vector<pathway> &PATHS,
                        const Rates           &rates,
                        const In_data         &input);

























#endif /* defined(__PumpKin__Print__) */
