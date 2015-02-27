//
//  pk_Print.cpp
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

#include "pk_Print.h"



// Printing individual pathway
void Print_pathway(const pathway &PATHWAY,
                   const In_data &input)
{
    if (PATHWAY.m_i.size() == 0)
    {
        cout << "There is no such a pathway! Choose other pathway to print." << endl;
        exit (EXIT_FAILURE);
    }
    cout << "---------------------------------------------" << endl;
    cout << "Reactions in the pathway:" << endl;
    for (size_t k = 0; k < PATHWAY.path.size(); k++)
    {
        if (PATHWAY.path[k] < n_R)
        {
            cout << PATHWAY.path[k] << ".  " << input.reactions[PATHWAY.path[k]] <<
                 "\t f = " << PATHWAY.f_k << ",\t x_j = " << PATHWAY.x_j[k] << endl;
        }
        else
        {
            // This is the species S of the pseudo-pathway ... -> S or S -> ...
            int s = PATHWAY.path[k] - n_R;

            if (PATHWAY.m_i[s] > 0)
            {
                cout << "     ... -> " << input.species[s] << endl;
            }
            else
            {
                cout << "     " << input.species[s] << " -> ..." << endl;
            }
        }
    }
    cout << endl;

    cout << "Number of molecules produced by pathway:" << endl;
    for (int i = 0; i < n_S; i++) if (PATHWAY.m_i[i] != 0)
            cout << input.species[i] << " -> " << PATHWAY.m_i[i] << endl;

    cout << "---------------------------------------------" << endl << endl;
}


// Printing all pathways
void Print_PATHS(const vector<pathway> &PATHS,
                 const Rates           &input)
{
    if (PATHS.size() == 0)
    {
        cout << "Nothing to print! PATHS is empty!" << endl;
        exit (EXIT_FAILURE);
    }

    doublearray1d print_p(PATHS.size());
    doublearray1d print_d(PATHS.size());

    print_p.setToValue(0);
    print_d.setToValue(0);

    for (int k = 0; k < PATHS.size(); k++)
    {
        if (PATHS[k].m_i[interest] > 0)
        {
            print_p(k) = print_p(k) + PATHS[k].m_i[interest] * PATHS[k].f_k;
            print_p(k) = 100.0 * print_p(k) / input.p_i(interest);
        }
        if (PATHS[k].m_i[interest] < 0)
        {
            print_d(k) = print_d(k) + PATHS[k].m_i[interest] * PATHS[k].f_k;
            print_d(k) = -100.0 * print_d(k) / input.d_i(interest);
        }
    }

    printf("+---------------------------------------------------+"
           "---------------------------+\n");
    printf("|       %-43s |       Rate of pathway     |\n",
           "Pathway");
    printf("+---------------------------------------------------+"
           "---------------------------+\n");

    for (size_t i = 0; i < PATHS.size(); i++)
    {
        cout.precision(3);
        if ((print_p(i) != 0.0) || (print_d(i) != 0.0) || true)
        {
            for (size_t j = 0; j < PATHS[i].path.size(); j++)
            {
                if (j == PATHS[i].path.size() / 2)
                {
                    printf("| %d * (%-36s) (%6d) | %-25g |\n",
                           PATHS[i].x_j[j],
                           input.reactions[PATHS[i].path[j]].c_str(), i,
                           PATHS[i].f_k);
                }
                else
                {
                    printf("| %d * (%-43s) |                           |\n",
                           PATHS[i].x_j[j],
                           input.reactions[PATHS[i].path[j]].c_str());
                }
            }
            printf("+---------------------------------------------------+"
                   "---------------------------+\n");
        }
    }

    printf("#####################################################"
           "############################\n");
    printf("+---------------------------------------------------+"
           "---------------------------+"
           "---------------------------+\n");
    printf("|       %-43s |  %-23s  |  %-24s |\n",
           "Procent of species from deleted pathways",
           "Production",
           "Consumption");
    printf("+---------------------------------------------------+"
           "---------------------------+"
           "---------------------------+\n");
    for (int i = 0; i < n_S; i++)
    {
        if ((input.tilda_p_i(i) * 100 / input.p_i(i) != 0.0) || (input.tilda_d_i(i) * 100 / input.d_i(i) != 0.0))
        {
            printf("| %-49s |"
                   "               %9.2g %% |"
                   "               %9.2g %% |\n",
                   input.species[i].c_str(),
                   (double) input.tilda_p_i(i) * 100 / input.p_i(i), (double) input.tilda_d_i(i) * 100 / input.d_i(i));

            printf("+---------------------------------------------------+"
                   "---------------------------+"
                   "---------------------------+\n");
        }
    }
    cout.precision(10);
}


// Printing double array in 2D
void Print_doublearray2D(const doublearray2d &A)
{
    for (int i = 0; i < A.getIndex1Size(); i++)
    {
        for (int j = 0; j < A.getIndex2Size(); j++)
            cout << A(i, j) << " ";
        cout << endl;
    }

}


// Printing final report
void Print_report(const vector<pathway> &PATHS,
                  const vector<int>     &all_brenching_points)
{
    cout << "Simulation ended successfully!" << endl << endl;
    cout << "Total number of pathways: " << PATHS.size() << endl;
    cout << "Total number of branching points: " << all_brenching_points.size() << endl;
}


// Printing all pathways according to specie of interest
void Print_PATHS_for_specie(const vector<pathway> &PATHS,
                            const Rates           &input)
{
    if (PATHS.size() == 0)
    {
        cout << "Nothing to print! PATHS is empty!" << endl;
        exit (EXIT_FAILURE);
    }

    doublearray1d print_p(PATHS.size());
    doublearray1d print_d(PATHS.size());

    print_p.setToValue(0);
    print_d.setToValue(0);

    for (int k = 0; k < PATHS.size(); k++)
    {
        if (PATHS[k].m_i[interest] > 0)
        {
            print_p(k) = print_p(k) + PATHS[k].m_i[interest] * PATHS[k].f_k;
            print_p(k) = 100.0 * print_p(k) / input.p_i(interest);
        }
        if (PATHS[k].m_i[interest] < 0)
        {
            print_d(k) = print_d(k) + PATHS[k].m_i[interest] * PATHS[k].f_k;
            print_d(k) = -100.0 * print_d(k) / input.d_i(interest);
        }
    }

    printf("+---------------------------------------------------+"
           "---------------------------+"
           "---------------------------+\n");
    printf("|       %-43s | Production of %-10s  | Consumption of %-10s |\n",
           "Pathway",
           input.species[interest].c_str(),
           input.species[interest].c_str());
    printf("+---------------------------------------------------+"
           "---------------------------+"
           "---------------------------+\n");

    for (size_t i = 0; i < PATHS.size(); i++)
    {
        cout.precision(3);
        if ((print_p(i) > 1.0) || (print_d(i) > 1.0) )
        {
            for (size_t j = 0; j < PATHS[i].path.size(); j++)
            {
                if (j == PATHS[i].path.size() / 2)
                {
                    printf("| %d * (%s) %22g |"
                           "               %9.2g %% |"
                           "               %9.2g %% |\n",
                           PATHS[i].x_j[j],
                           input.reactions[PATHS[i].path[j]].c_str(),
                           PATHS[i].f_k,
                           (double) print_p(i), (double) print_d(i));
                }
                else
                {
                    printf("| %d * (%-35s)         |"
                           "                           |"
                           "                           |\n",
                           PATHS[i].x_j[j],
                           input.reactions[PATHS[i].path[j]].c_str());
                }
            }
            printf("+---------------------------------------------------+"
                   "---------------------------+"
                   "---------------------------+\n");
        }
    }

    printf("#####################################################"
           "############################"
           "############################\n");
    printf("+---------------------------------------------------+"
           "---------------------------+"
           "---------------------------+\n");
    printf("|       %-43s |  %-23s  |  %-24s |\n",
           "Procent of species from deleted pathways",
           "Production",
           "Consumption");
    printf("+---------------------------------------------------+"
           "---------------------------+"
           "---------------------------+\n");
    for (int i = 0; i < n_S; i++)
    {
        if ((input.tilda_p_i(i) * 100 / input.p_i(i) != 0.0) || (input.tilda_d_i(i) * 100 / input.d_i(i) != 0.0))
        {
            printf("| %-49s |"
                   "               %9.2g %% |"
                   "               %9.2g %% |\n",
                   input.species[i].c_str(),
                   (double) input.tilda_p_i(i) * 100 / input.p_i(i), (double) input.tilda_d_i(i) * 100 / input.d_i(i));

            printf("+---------------------------------------------------+"
                   "---------------------------+"
                   "---------------------------+\n");
        }
    }
    cout.precision(10);

}


// Printing all reduced system
void Print_PATHS_reduce(const vector<pathway> &PATHS,
                        const Rates           &rates,
                        const In_data         &input)

{
    if (PATHS.size() == 0)
    {
        cout << "Nothing to print! PATHS is empty!" << endl;
        exit (EXIT_FAILURE);
    }

    cout << "Aram \n" << endl;
    for (int i = 0; i < PATHS.size(); i++)
    {
        cout << "Pathway number " << i << endl;
        Print_pathway(PATHS[i], input);
    }
    for (int j = 0; j < PATHS.size(); j++)
    {
        for (int i = 0; i < input.species.size(); i++)
        {
            cout << PATHS[j].m_i[i]*PATHS[j].f_k / rates.d_i(i) << "     ";
        }
        cout << endl;
    }

}

























