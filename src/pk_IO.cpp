//
//  pk_IO.cpp
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

#include "pk_IO.h"

// Print license information
void Print_license()
{
    cout << "+++**********************************************************************+++" << endl;
    cout << "+++                                                                      +++" << endl;
    cout << "+++ PumpKin: A tool to find principal pathways in plasma chemical models +++" << endl;
    cout << "+++ Copyright\u00A9 2013-2015 Aram H Markosyan.                               +++" << endl;
    cout << "+++                                                                      +++" << endl;
    cout << "+++ This program is free software; you can redistribute it and/or        +++" << endl;
    cout << "+++ modify it under the terms of the GNU General Public License          +++" << endl;
    cout << "+++ as published by the Free Software Foundation; either version 2       +++" << endl;
    cout << "+++ of the License, or (at your option) any later version.               +++" << endl;
    cout << "+++ This program is distributed in the hope that it will be useful,      +++" << endl;
    cout << "+++ but WITHOUT ANY WARRANTY; without even the implied warranty of       +++" << endl;
    cout << "+++ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +++" << endl;
    cout << "+++ GNU General Public License for more details.                         +++" << endl;
    cout << "+++                                                                      +++" << endl;
    cout << "+++ Point of Contact:   Dr. Aram H. Markosyan                            +++" << endl;
    cout << "+++ Address: University of Michigan, Electrical Engineering and Computer +++" << endl;
    cout << "+++ Science Department, 1301 Beal Ave, Ann Arbor, MI 48109-2122          +++" << endl;
    cout << "+++ Email: armarkos@umich.edu                                            +++" << endl;
    cout << "+++ Tel: 734-647-4840                                                    +++" << endl;
    cout << "+++ Homepage: http://markosyanaram.com                                   +++" << endl;
    cout << "+++                                                                      +++" << endl;
    cout << "+++**********************************************************************+++" << endl;

    cout << endl << endl;
}

// Reads input kinetic model
void Read_kin(In_data      &kinetics,
              const string &folder)
{
    string m_folder;

    if (folder == "")
    {
        ifstream input_file("Input/input.txt");
        if (!input_file)
        {
            input_file.close();
            input_file.open("/Users/Apple/Devel/PumpKin/Input/input.txt");
            if (input_file)
            {
                m_folder = "/Users/Apple/Devel/PumpKin/Input/";
                input_file.close();
            }
            else
            {
                /* could not open directory */
                perror ("Input folder or input file");
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            input_file.close();
            m_folder = "Input/";
        }
    }
    else
    {
        if (folder.at(folder.size() - 1) == '/') m_folder = folder;
        else m_folder = folder + "/";
    }

    string filename = m_folder + "input.txt";
    ifstream file(filename.c_str());

    if (!file)
    {
        // file couldn't be opened
        cerr << "Error: File containing the parameters could not be opened" << endl;
        exit(1);
    }

    string line;
    for (int i = 0; i < 8; i++)
    {
        getline(file, line);
        string temp1;
        istringstream iss(line);
        iss >> temp1;
        iss >> temp1;
        if (i == 0) iss >> interest;
        if (i == 1) iss >> t_init;
        if (i == 2) iss >> t_end;
        if (i == 3) iss >> max_bp;
        if (i == 4) iss >> tau_lifetime;
        if (i == 5) iss >> max_path;
        if (i == 6) iss >> f_min;
        if (i == 7) iss >> global_kin;
    }

    file.close();
    interest = interest - 1; // to align with C++ notations

    if (t_end - t_init < 0)
    {
        cout << "Error: Wrong time period: t_end < t_init." << endl;
        exit (EXIT_FAILURE);
    }

    Read_species(kinetics, m_folder);
    Read_reactions(kinetics, m_folder);
    Read_time(kinetics, m_folder);
    Read_density(kinetics, m_folder);
    Read_rates(kinetics, m_folder);
    Read_matrix(kinetics, m_folder);
}

// Reads species from files
void Read_species(In_data      &kinetics,
                  const string &m_folder)
{
    if (global_kin == false)
    {
        ifstream file;
        string filename_1 = m_folder + "qt_species_list.txt";
        string filename_2 = m_folder + "QT_SPECIES_LIST.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: File containing the list of species could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());

        // Read file line by line
        string line;
        while (getline(file, line))
        {
            string temp = line.substr(line.find(" ", 1) + 1);
            temp.erase(std::find_if(temp.rbegin(), temp.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), temp.end());

            if (temp != "")
                kinetics.species.push_back(temp);
        }

        n_S = kinetics.species.size();
        file.close();
    }
    else
    {
        ifstream file;
        string filename_1 = m_folder + "pumpkin_species_list.txt";
        string filename_2 = m_folder + "PUMPKIN_SPECIES_LIST.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: File containing the list of species could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());

        // Read file line by line
        string line;
        while (getline(file, line))
        {
            string temp = line.substr(line.find(" ", 1) + 1);
            if ( line.size() < 2 )
                break;

            // Removing trailing spaces
            temp.erase(find_if(temp.rbegin(), temp.rend(), bind1st(not_equal_to<char>(), ' ')).base(), temp.end());

            if (temp != "")
                kinetics.species.push_back(temp);
        }

        n_S = kinetics.species.size();
        file.close();
    }

    if ((interest > 0) && (n_S < interest))
    {
        cerr << "Error: The species of interet is not in your chemistry file." << endl;
        cerr << "You have only " << n_S << " species." << endl;
        exit(1);
    }
}

// Reads the stoichiometric matrix from files
void Read_matrix(In_data      &kinetics,
                 const string &m_folder)
{
    if (global_kin == false)
    {
        ifstream file;
        string filename_1 = m_folder + "qt_matrix.txt";
        string filename_2 = m_folder + "QT_MATRIX.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: File containing the stoichiometric matrix could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());

        string line;
        int i = 0;
        n_R = 100000000;
        while (getline(file, line))
        {
            if (!Check_file_line(line))
            {
                cout << "Error: The stoichiometric matrix contains inf values!" << endl;
                exit (EXIT_FAILURE);
            }

            int j = 0;
            istringstream iss(line);
            int val;

            if (i == 0)
            {
                vector<int> temp;
                while (iss >> val)
                {
                    temp.push_back(val);
                    j++;
                }
                n_R = j;
                kinetics.matrix.initialize(n_S, n_R);
                for (int k = 0; k < n_R; k++) kinetics.matrix(i, k) = temp[k];
            }
            else
            {
                while (iss >> val)
                {
                    kinetics.matrix(i, j) = val;
                    j++;
                }
            }
            i++;
            if (n_R == 100000000) n_R = j;
            else
            {
                if (n_R != j)
                {
                    cerr << "Something is wrong with stoichiometric matrix" << endl;
                    exit(1);
                }
            }
        }


        file.close();
    }
    else
    {
        ifstream file;
        string filename_1 = m_folder + "pumpkin_matrix.txt";
        string filename_2 = m_folder + "PUMPKIN_MATRIX.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: File containing the stoichiometric matrix could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());

        string line;
        double temp;
        kinetics.matrix.initialize(n_S, n_R);
        for (int i = 0; i < n_S; i++)
        {
            for (int j = 0; j < n_R; j++)
            {
                getline(file, line);
                if (!Check_file_line(line))
                {
                    cout << "Error: The the file containing densities contains non-integer values!" << endl;
                    exit (EXIT_FAILURE);
                }
                else if ( line.size() < 2 ) break;
                istringstream iss(line);
                iss >> kinetics.matrix(i, j);
            }
        }


        file.close();
    }
}


// Checking if all values are numbers
bool Check_file_line(string &line)
{
    size_t nan_pos = line.find("nan");

    if (nan_pos != string::npos)
    {
        cout << "Warning: Densities matrix contain nan value. It is replaced by 0." << endl;
        line.replace(nan_pos, 3, "0.0");
        return true;
    }
    else if (line.find("inf") != string::npos)
        return false;
    return true;
}


// Reads the time steps from files
void Read_time(In_data      &kinetics,
               const string &m_folder)
{
    if (global_kin == false)
    {
        ifstream file;
        string filename_1 = m_folder + "qt_conditions.txt";
        string filename_2 = m_folder + "QT_CONDITIONS.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: File containing time could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());

        file.ignore(2000, '\n');

        string line;
        int i = 0;
        vector<double> temp;
        while (getline(file, line))
        {
            if (!Check_file_line(line))
            {
                cout << "Error: The the file containing time contains inf values!" << endl;
                exit (EXIT_FAILURE);
            }
            else if ( line.size() < 2 ) break;
            istringstream iss(line);
            double val;
            iss >> val;
            temp.push_back(val);
            i++;
        }

        n_t = i;
        kinetics.time.initialize(n_t);
        for (int i = 0; i < n_t; i++) kinetics.time(i) = temp[i];
        file.close();
    }
    else
    {
        ifstream file;
        string filename_1 = m_folder + "pumpkin_conditions.txt";
        string filename_2 = m_folder + "PUMPKIN_CONDITIONS.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: File containing time could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());

        string line;
        int i = 0;
        vector<double> temp;
        while (getline(file, line))
        {
            if (!Check_file_line(line))
            {
                cout << "Error: The the file containing time contains inf values!" << endl;
                exit (EXIT_FAILURE);
            }
            else if ( line.size() < 2 ) break;
            istringstream iss(line);
            double val;
            iss >> val;
            temp.push_back(val);
            i++;
        }

        n_t = i;
        kinetics.time.initialize(n_t);
        for (int i = 0; i < n_t; i++) kinetics.time(i) = temp[i];
        file.close();
    }

}


// Reads the density from files
void Read_density(In_data      &kinetics,
                  const string &m_folder)
{
    if (global_kin == false)
    {
        ifstream file;
        string filename_1 = m_folder + "qt_densities.txt";
        string filename_2 = m_folder + "QT_DENSITIES.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: File containing densities could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());

        kinetics.densities.resize(n_t);

        file.ignore(200000, '\n');

        string line;
        double temp;
        for (int i = 0; i < n_t; i++)
        {
            kinetics.densities[i].resize(n_S);
            getline(file, line);
            if (!Check_file_line(line))
            {
                cout << "Error: The the file containing densities contains non-integer values!" << endl;
                exit (EXIT_FAILURE);
            }
            istringstream iss(line);
            iss >> temp;
            int j = 0;
            while (j < n_S)
            {
                iss >> kinetics.densities[i](j);
                j++;
            }
        }

        file.close();
    }
    else
    {
        ifstream file;
        string filename_1 = m_folder + "pumpkin_densities.txt";
        string filename_2 = m_folder + "PUMPKIN_DENSITIES.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: File containing densities could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());

        kinetics.densities.resize(n_t);

        string line;
        double temp;
        for (int i = 0; i < n_t; i++)
        {
            kinetics.densities[i].resize(n_S);
            for (int j = 0; j < n_S; j++)
            {
                getline(file, line);
                if (!Check_file_line(line))
                {
                    cout << "Error: The the file containing densities contains non-integer values!" << endl;
                    exit (EXIT_FAILURE);
                }
                else if ( line.size() < 2 ) break;
                istringstream iss(line);
                iss >> kinetics.densities[i](j);
            }
        }

        file.close();
    }

}


// Reads the the rates from files
void Read_rates(In_data      &kinetics,
                const string &m_folder)
{
    if (global_kin == false)
    {
        ifstream file;
        string filename_1 = m_folder + "qt_rates.txt";
        string filename_2 = m_folder + "QT_RATES.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: The file containing rates could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());

        kinetics.rates.resize(n_t);

        file.ignore(200000, '\n');
        string line;
        double temp;
        for (int i = 0; i < n_t; i++)
        {
            kinetics.rates[i].resize(n_R);
            getline(file, line);
            if (!Check_file_line(line))
            {
                cout << "Error: The file containing rates contains non-integer values!" << endl;
                exit (EXIT_FAILURE);
            }
            istringstream iss(line);
            iss >> temp;
            int j = 0;
            while (j < n_R)
            {
                iss >> kinetics.rates[i](j);
                j++;
            }
        }


        file.close();
    }
    else
    {
        ifstream file;
        string filename_1 = m_folder + "pumpkin_rates.txt";
        string filename_2 = m_folder + "PUMPKIN_RATES.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: The file containing rates could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());

        kinetics.rates.resize(n_t);

        string line;
        double temp;
        for (int i = 0; i < n_t; i++)
        {
            kinetics.rates[i].resize(n_R);
            for (int j = 0; j < n_R; j++)
            {
                getline(file, line);
                if (!Check_file_line(line))
                {
                    cout << "Error: The the file containing densities contains non-integer values!" << endl;
                    exit (EXIT_FAILURE);
                }
                else if ( line.size() < 2 ) break;
                istringstream iss(line);
                iss >> kinetics.rates[i](j);
            }
        }

        file.close();
    }

}


// Reads the the reactions from files
void Read_reactions(In_data      &kinetics,
                    const string &m_folder)
{
    if (global_kin == false)
    {
        ifstream file;
        string filename_1 = m_folder + "qt_reactions_list.txt";
        string filename_2 = m_folder + "QT_REACTIONS_LIST.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: File containing the list of reactions could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());


        int i = 1;
        string line, temp;
        while (getline(file, line))
        {
            stringstream ss; ss << i; string index = ss.str();
            temp = line.substr(line.find(index) + index.size() + 1);
            //string temp = line.substr(line.find(" ")+1);
            temp.erase(std::find_if(temp.rbegin(), temp.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), temp.end());


            if (temp != "") kinetics.reactions.push_back(temp);
            i++;

            temp = "";
        }

        n_R = kinetics.reactions.size();
        file.close();
    }
    else
    {
        ifstream file;
        string filename_1 = m_folder + "pumpkin_reaction_list.txt";
        string filename_2 = m_folder + "PUMPKIN_REACTION_LIST.TXT"; // In case of VAX/VMS Operating systems
        ifstream file_1(filename_1.c_str());
        if (!file_1)
        {
            ifstream file_2(filename_2.c_str());
            if (!file_2)
            {
                // file couldn't be opened
                cerr << "Error: File containing the list of reactions could not be opened" << endl;
                exit(1);
            }
            else file.open(filename_2.c_str());
        }
        else file.open(filename_1.c_str());


        int i = 1;
        string line, temp;
        while (getline(file, line))
        {
            temp = line.substr(0);
            if ( line.find(">") == string::npos ) break;

            // Removing trailing spaces
            temp.erase(find_if(temp.rbegin(), temp.rend(), bind1st(not_equal_to<char>(), ' ')).base(), temp.end());

            if (temp != "") kinetics.reactions.push_back(temp);
            i++;

            temp = "";
        }

        n_R = kinetics.reactions.size();
        file.close();
    }

}

// Average input kinetic model
void Average_kin(In_data &kinetics,
                 Rates   &rates)
{
    double delta_t = t_end - t_init;

    if (delta_t < 0)
    {
        cout << "Error: Wrong time period: t_end < t_init." << endl;
        exit (EXIT_FAILURE);
    }

    int i_init = 0, i_end = 0;

    for (int i = 0; i < n_t; i++)
    {
        i_init = i;
        if (kinetics.time(i) >= t_init) break;
    }

    for (int i = 0; i < n_t - 1; i++)
    {
        i_end = i;
        if (kinetics.time(i) >= t_end) break;

    }

    doublearray1d r_init(n_R), r_end(n_R);
    interp_1(r_init, t_init, kinetics.time, kinetics.rates);
    interp_1(r_end, t_end, kinetics.time, kinetics.rates);

    doublearray1d c_init(n_S), c_end(n_S);
    interp_1(c_init, t_init, kinetics.time, kinetics.densities);
    interp_1(c_end, t_end, kinetics.time, kinetics.densities);

    rates.r_j.initialize(n_R);
    for (int j = 0; j < n_R; j++) rates.r_j(j) = 0.0;
    if (t_init == t_end) for (int j = 0; j < n_R; j++)
            rates.r_j(j) = r_init(j);
    else
    {
        for (int j = 0; j < n_R; j++)
        {
            rates.r_j(j) = (kinetics.time(i_init + 1) - t_init) * 0.5 * (r_init(j) + kinetics.rates[i_init + 1](j));

            for (int i = i_init + 1; i < i_end - 1; i++)
                rates.r_j(j) = rates.r_j(j) + (kinetics.time(i + 1) - kinetics.time(i)) * 0.5 * (kinetics.rates[i](j) + kinetics.rates[i + 1](j));

            if (i_init != i_end - 1)
            {
                rates.r_j(j) = rates.r_j(j) + (t_end - kinetics.time(i_end - 1)) * 0.5 * (r_end(j) + kinetics.rates[i_end - 1](j));
            }

            rates.r_j(j) = rates.r_j(j) / delta_t;
        }

    }

    rates.c_i.initialize(n_S);
    for (int i = 0; i < n_S; i++) rates.c_i(i) = 0.0;
    if (t_init == t_end) for (int j = 0; j < n_S; j++)
            rates.c_i(j) = c_init(j);
    else
    {
        for (int j = 0; j < n_S; j++)
        {
            rates.c_i(j) = (kinetics.time(i_init + 1) - t_init) * 0.5 * (c_init(j) + kinetics.densities[i_init + 1](j));

            for (int i = i_init + 1; i < i_end - 1; i++)
                rates.c_i(j) = rates.c_i(j) + (kinetics.time(i + 1) - kinetics.time(i)) * 0.5 * (kinetics.densities[i](j) + kinetics.densities[i + 1](j));

            if (i_init != i_end - 1)
            {
                rates.c_i(j) = rates.c_i(j) + (t_end - kinetics.time(i_end - 1)) * 0.5 * (c_end(j) + kinetics.densities[i_end - 1](j));
            }

            rates.c_i(j) = rates.c_i(j) / delta_t;
        }

    }

    bool temp = true;
    rates.delta_c.initialize(n_S);
    for (int i = 0; i < n_S; i++)
    {
        rates.delta_c(i) = 0;
        for (int j = 0; j < n_R; j++)
            rates.delta_c(i) = rates.delta_c(i) + kinetics.matrix(i, j) * rates.r_j(j) * delta_t;

        //temp = (((c_end(i) - c_init(i)) - rates.delta_c(i)) == 0);
    }

    // if (!temp) cout<<"Worrning: the input files violate conservation!\n"<<endl;

    for (int i = 0; i < n_S; i++) rates.species.push_back(kinetics.species[i]);
    for (int j = 0; j < n_R; j++) rates.reactions.push_back(kinetics.reactions[j]);

}


// Interpolation function
void interp_1(doublearray1d               &interpol,
              const double                &t,
              const doublearray1d         &T,
              const vector<doublearray1d> &Y)
{
    int j = 0;

    if (t < T(0))
    {
        cout << "Error: Wrong time period: t_init is out of kinetic model's time bounds." << endl;
        cout << "Time zero is: " << T(0) << "\n";
        exit (EXIT_FAILURE);
    }

    if (t > T(n_t - 1))
    {
        cout << "Error: Wrong time period: t_end is out of kinetic model's time bounds." << endl;
        cout << "End time is: " << T(n_t - 1) << "\n";
        exit (EXIT_FAILURE);
    }

    while (j <= n_t - 1)
    {
        if (T(j) >= t) break;
        j = j + 1;
    }

    if (t == T(0)) for (int i = 0; i < Y[0].getIndex1Size(); i++)
            interpol(i) = Y[0](i);
    else if (t == T(n_t - 1)) for (int i = 0; i < Y[n_t - 1].getIndex1Size(); i++)
            interpol(i) = Y[n_t - 1](i);
    else for (int i = 0; i < Y[0].getIndex1Size(); i++)
            interpol(i) = Y[j - 1](i) + (Y[j](i) - Y[j - 1](i)) * (t - T(j - 1)) / (T(j) - T(j - 1));
}































