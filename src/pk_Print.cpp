//
//  pk_Print.cpp
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

// Convert Global_Kin results to ZDPlasKin results
void Global_Kin__to_ZDPlasKin(const In_data &kinetics,
                              string        &folder)
{
  folder = folder + "/QtPlasKin";

  ofstream file;
  string filename;

  // Creating the file "qt_conditions_list.txt"
  filename = folder + "/qt_conditions_list.txt";
  file.open(filename.c_str());

  if (!file.is_open()) {
    cout << "Can't open the file: qt_conditions_list.txt\n";
    cout << "Create directory: " << folder << "\n";
    exit(1);
  }
  else {
    file << "  1 Reduced field [Td]\n";
    file << "  2 Gas temperature [K]\n";
    file << "  3 Electron temperature [K]\n";
    file << "  4 Current density [A/cm2]\n";
    file << "  5 Power density [W/cm3]\n";
  }

  file.close();

  // Creating the file "qt_species_list.txt"
  filename = folder + "/qt_species_list.txt";
  file.open(filename.c_str());

  if (!file.is_open()) {
    cout << "Can't open the file: qt_species_list.txt\n";
    cout << "Create directory: " << folder << "\n";
    exit(1);
  }
  else {
    for (size_t i = 0; i < kinetics.species.size(); i++) {
      if (kinetics.species.size() > 100) {
        if (i < 10)
          file << "  " << i + 1 << " " << kinetics.species[i] << "\n";
        else if (i < 100)
          file << " " << i + 1 << " " << kinetics.species[i] << "\n";
        else
          file << i + 1 << " " << kinetics.species[i] << "\n";
      } else {
        if (i < 10)
          file << " " << i + 1 << " " << kinetics.species[i] << "\n";
        else
          file << i + 1 << " " << kinetics.species[i] << "\n";
      }
    }
  }

  file.close();

  // Creating the file "qt_reactions_list.txt"
  filename = folder + "/qt_reactions_list.txt";
  file.open(filename.c_str());

  if (!file.is_open()) {
    cout << "Can't open the file: qt_reactions_list.txt\n";
    cout << "Create directory: " << folder << "\n";
    exit(1);
  }
  else {
    for (size_t i = 0; i < kinetics.reactions.size(); i++) {
      if (kinetics.reactions.size() > 1000) {
        if (i < 10)
          file << "   " << i + 1 << " " << kinetics.reactions[i] << "\n";
        if (i < 100)
          file << "  " << i + 1 << " " << kinetics.reactions[i] << "\n";
        if (i < 1000)
          file << " " << i + 1 << " " << kinetics.reactions[i] << "\n";
        else
          file << i + 1 << " " << kinetics.reactions[i] << "\n";
      } else if (kinetics.reactions.size() > 100) {
        if (i < 10)
          file << "  " << i + 1 << " " << kinetics.reactions[i] << "\n";
        else if (i < 100)
          file << " " << i + 1 << " " << kinetics.reactions[i] << "\n";
        else
          file << i + 1 << " " << kinetics.reactions[i] << "\n";
      } else {
        if (i < 10)
          file << " " << i + 1 << " " << kinetics.reactions[i] << "\n";
        else
          file << i + 1 << " " << kinetics.reactions[i] << "\n";
      }
    }
  }

  file.close();

  // Creating the file "qt_matrix.txt"
  filename = folder + "/qt_matrix.txt";
  file.open(filename.c_str());

  if (!file.is_open()) {
    cout << "Can't open the file: qt_matrix.txt\n";
    cout << "Create directory: " << folder << "\n";
    exit(1);
  }
  else {
    for (int i = 0; i < kinetics.species.size(); i++) {
      for (int j = 0; j < kinetics.reactions.size(); j++) {
        file << kinetics.matrix(i, j) << "\t";
      }
      file << "\n";
    }
  }

  file.close();

  // Creating the file "qt_densities.txt"
  filename = folder + "/qt_densities.txt";
  file.open(filename.c_str());

  if (!file.is_open()) {
    cout << "Can't open the file: qt_densities.txt\n";
    cout << "Create directory: " << folder << "\n";
    exit(1);
  }
  else {
    file << "   Time_s";
    for (int i = 0; i < kinetics.species.size(); i++)
      file << "\t" << i + 1;
    file << "\n";

    for (int i = 0; i < kinetics.time.getIndex1Size(); i++) {
      file.setf(ios::scientific | ios::showpoint);
      cout.precision(6);
      file << kinetics.time(i);
      for (int j = 0; j < kinetics.species.size(); j++) {
        file << "\t" << kinetics.densities[i](j);
      }
      file << "\n";
    }
  }

  file.close();

  // Creating the file "qt_rates.txt"
  filename = folder + "/qt_rates.txt";
  file.open(filename.c_str());

  if (!file.is_open()) {
    cout << "Can't open the file: qt_rates.txt\n";
    cout << "Create directory: " << folder << "\n";
    exit(1);
  }
  else {
    file << "   Time_s";
    for (int i = 0; i < kinetics.reactions.size(); i++)
      file << "\t" << i + 1;
    file << "\n";

    for (int i = 0; i < kinetics.time.getIndex1Size(); i++) {
      file.setf(ios::scientific | ios::showpoint);
      cout.precision(6);
      file << kinetics.time(i);
      for (int j = 0; j < kinetics.reactions.size(); j++) {
        file << "\t" << kinetics.rates[i](j);
      }
      file << "\n";
    }
  }

  file.close();

  // Creating the file "qt_conditions.txt"
  filename = folder + "/qt_conditions.txt";
  file.open(filename.c_str());

  if (!file.is_open()) {
    cout << "Can't open the file: qt_conditions.txt\n";
    cout << "Create directory: " << folder << "\n";
    exit(1);
  }
  else {
    file << "   Time_s";
    for (int i = 0; i < 5; i++)
      file << "\t" << i + 1;
    file << "\n";

    for (int i = 0; i < kinetics.time.getIndex1Size(); i++) {
      file.setf(ios::scientific | ios::showpoint);
      cout.precision(6);
      file << kinetics.time(i);
      file << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0;
      file << "\n";
    }
  }

  file.close();
}

// Printing all pathways
void Print_PATHS(const vector<pathway> &PATHS,
                 const Rates           &input,
                 const int             &interst_species)
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
    if (PATHS[k].m_i[interst_species] > 0)
    {
      print_p(k) = print_p(k) + PATHS[k].m_i[interst_species] * PATHS[k].f_k;
      print_p(k) = 100.0 * print_p(k) / input.p_i(interst_species);
    }
    if (PATHS[k].m_i[interst_species] < 0)
    {
      print_d(k) = print_d(k) + PATHS[k].m_i[interst_species] * PATHS[k].f_k;
      print_d(k) = -100.0 * print_d(k) / input.d_i(interst_species);
    }
  }

  printf("+---------------------------------------------------+"
         "---------------------------+\n");
  printf("|       %-43s |       Rate of pathway     |\n",
         "Pathway");
  printf("+---------------------------------------------------+"
         "---------------------------+\n");

  for (int i = 0; i < PATHS.size(); i++)
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
  cout << "\n" << "Simulation ended successfully!" << endl << endl;
  cout << "Total number of pathways: " << PATHS.size() << endl;
  cout << "Total number of branching points: " << all_brenching_points.size() << endl;
}


// Printing all pathways according to specie of interest
void Print_PATHS_for_specie(const vector<pathway> &PATHS,
                            const Rates           &input,
                            const int             &interst_species)
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
    if (PATHS[k].m_i[interst_species] > 0)
    {
      print_p(k) = print_p(k) + PATHS[k].m_i[interst_species] * PATHS[k].f_k;
      print_p(k) = 100.0 * print_p(k) / input.p_i(interst_species);
    }
    if (PATHS[k].m_i[interst_species] < 0)
    {
      print_d(k) = print_d(k) + PATHS[k].m_i[interst_species] * PATHS[k].f_k;
      print_d(k) = -100.0 * print_d(k) / input.d_i(interst_species);
    }
  }

  printf("+---------------------------------------------------+"
         "---------------------------+"
         "---------------------------+\n");
  printf("|       %-43s | Production of %-10s  | Consumption of %-10s |\n",
         "Pathway",
         input.species[interst_species].c_str(),
         input.species[interst_species].c_str());
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


// Print interactive reports
void Print_Results(const vector<pathway> &PATHS,
                   const Rates           &rates,
                   const In_data         &kinetics,
                   string                &m_folder,
                   const vector<int>     &all_brenching_points)
{
  if (interest < 0)
    Print_PATHS(PATHS, rates, interest);
  else
    Print_PATHS_for_specie(PATHS, rates, interest);

  bool more_interest = true;
  int new_interest;

  while (more_interest) {
    cout << "\n";
    printf("#####################################################"
           "############################"
           "#####################################\n");
    cout << "\n" << "Enter positive for a species of interest, 0 for dominant pathways, negative number for termination: ";
    cin >> new_interest;

    if (new_interest < 0) {
      if (new_interest == -1000) Global_Kin__to_ZDPlasKin(kinetics, m_folder);
      more_interest = false;
    }
    else if (new_interest > rates.species.size()) {
      cout << "\n" << "Maximum number of species in your system is: " << rates.species.size() << "\n";
      cout << "\n" << "Enter positive for a species of interest, 0 for dominant pathways, negative number for termination: ";
      cin >> new_interest;
      if (new_interest < 0)
        more_interest = false;
      else if (new_interest > rates.species.size())
        more_interest = false;
      else if (new_interest > 0) {
        Print_PATHS_for_specie(PATHS, rates, new_interest);
        Print_report(PATHS, all_brenching_points);
        more_interest = true;
      }
      else if (new_interest == 0) {
        Print_PATHS(PATHS, rates, new_interest);
        Print_report(PATHS, all_brenching_points);
        more_interest = true;
      }
      else
        more_interest = false;
    }
    else if (new_interest > 0) {
      Print_PATHS_for_specie(PATHS, rates, new_interest);
      Print_report(PATHS, all_brenching_points);
      more_interest = true;
    }
    else if (new_interest == 0) {
      Print_PATHS(PATHS, rates, new_interest);
      Print_report(PATHS, all_brenching_points);
      more_interest = true;
    }
    else
      more_interest = false;
  }

  cout << "\n" << "Good bye!\n\n";
}






















