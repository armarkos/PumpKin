//
//  pk_Pathways.cpp
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

#include "pk_Pathways.h"

// Initialization of pathways
void Initialize_PATHS(vector<pathway> &PATHS,
                      const In_data   &kinetics,
                      const Rates     &rates)
{
    for (int j = 0; j < n_R; j++)
    {
        PATHS[j].deleted = false;
        PATHS[j].simple = true;
        PATHS[j].path.push_back(j);
        PATHS[j].x_j.push_back(1);
        PATHS[j].f_k = rates.r_j(j);
        for (int i = 0; i < n_S; i++)
            PATHS[j].m_i.push_back(kinetics.matrix(i, j));
    }
}


// Calculate rates
void Initialize_rates(vector<pathway> &PATHS,
                      const In_data   &kinetics,
                      Rates           &rates)
{
    rates.p_i.initialize(n_S);
    rates.d_i.initialize(n_S);
    rates.D_i.initialize(n_S);
    rates.delta_i.initialize(n_S);
    rates.d_i.initialize(n_S);

    rates.tilda_d_i.initialize(n_S);
    rates.tilda_p_i.initialize(n_S);
    rates.tilda_r_j.initialize(n_R);

    for (int i = 0; i < n_S; i++)
    {
        rates.tilda_d_i(i) = 0.0;
        rates.tilda_p_i(i) = 0.0;
        rates.tilda_r_j(i) = 0.0;
    }

    Calculate_rates(PATHS, kinetics, rates);

}


// Find and sort brenching points
void Branching_points(vector<int> &bp,
                      const Rates &rates)
{
    for (int i = 0; i < n_S; i++)
        if (isfinite(rates.c_i(i) / rates.d_i(i)))
            if (rates.c_i(i) / rates.d_i(i) > tau_lifetime)
            {
                bp.push_back(i);
            }

    // Cheking whether tau_lifetime is not too short
    if (rates.c_i(interest) / rates.d_i(interest) < tau_lifetime)
    {
        cout << "The lifetime of the specie " << rates.species[interest] << " is shorter than the 'minimal life time'!" << endl << "Decrease tau_lifetime to less than the value " << rates.c_i(interest) / rates.d_i(interest) << endl;
        exit (EXIT_FAILURE);
    }

    // Sorting all_brenching_points
    for (int i = 0; i < bp.size(); i++)
    {
        for (int j = i + 1; j < bp.size(); j++)
        {
            if (rates.c_i(bp[i]) / rates.d_i(bp[i]) > rates.c_i(bp[j]) / rates.d_i(bp[j]))
            {
                int t = bp[i];
                bp[i] = bp[j];
                bp[j] = t;
            }
        }
    }
}

// Get the next brenching point, returns false if there is no available one
bool Get_bp(vector<int> &all_bp,
            const Rates &rates,
            int         &bp)
{
    double temp = 1.0e40;
    bp = -5;

    if (max_bp <= 0)
    {
        for (int i = 0; i < n_S; i++) //if (i != interest)
            if (isfinite(rates.c_i(i) / rates.d_i(i))) // is finite
                if ((rates.c_i(i) / rates.d_i(i) <= tau_lifetime) && (rates.c_i(i) / rates.d_i(i) < temp)) // smaller than lifetime and current one
                    if (std::find(all_bp.begin(), all_bp.end(), i) == all_bp.end())   // has not been brenching point so far
                    {
                        temp = rates.c_i(i) / rates.d_i(i);
                        bp = i;
                    }
        if (bp == -5) return false;
        else
        {
            all_bp.push_back(bp);
            return true;
        }
    }
    else
    {
        if (all_bp.size() < max_bp)
        {
            for (int i = 0; i < n_S; i++) //if (i != interest)
                if (isfinite(rates.c_i(i) / rates.d_i(i))) // is finite
                    if ((rates.c_i(i) / rates.d_i(i) <= tau_lifetime) && (rates.c_i(i) / rates.d_i(i) < temp)) // smaller than lifetime and current one
                        if (std::find(all_bp.begin(), all_bp.end(), i) == all_bp.end())   // has not been brenching point so far
                        {
                            temp = rates.c_i(i) / rates.d_i(i);
                            bp = i;
                        }
            if (bp == -5) return false;
            else
            {
                //cout << "Lifetime of " << bp << " " << temp << " s" << endl;
                all_bp.push_back(bp);
                return true;
            }
        }
        else
            return false;
    }
}


// Run baby, run!
void Run(vector<pathway>   &PATHS,
         vector<int>       &bp,
         const In_data     &kinetics,
         Rates             &rates)
{
    int brenching_point;
    while (Get_bp(bp, rates, brenching_point))
    {
        cout << bp.size() << ". Now treating branching point " << rates.species[brenching_point] << ", with lifetime " << rates.c_i(brenching_point) / rates.d_i(brenching_point) << ", is treating. The number of pathways is " << PATHS.size() << "." << endl;
        Merge_pathways(PATHS, brenching_point, kinetics, rates);
        Insignificant(PATHS, rates, kinetics);

        // Recalculate the p, d
        Calculate_rates(PATHS, kinetics, rates);

        int count = 0;
        for (int i = 0; i < PATHS.size(); i++)
        {
            if (!PATHS[i].deleted) count++;
        }

        vector<elementary> Elementary(PATHS.size());
        Make_elementary(PATHS, Elementary, bp, rates, kinetics);

        Splitting(PATHS, Elementary, bp, rates, kinetics);
    }
}


// Merge pathways at this brenching point
void Merge_pathways(vector<pathway> &PATHS,
                    const int       &brenching_point,
                    const In_data   &kinetics,
                    Rates           &rates)
{
    size_t size_temp = PATHS.size();
    double dtilde = rates.tilda_d_i(brenching_point);
    double ptilde = rates.tilda_p_i(brenching_point);

    vector<int> PATHS_producing;
    for (size_t i = 0; i < size_temp; i++)
    {
        if (PATHS[i].m_i[brenching_point] > 0)
            PATHS_producing.push_back((int)i);
    }

    vector<int> PATHS_consuming;
    for (size_t i = 0; i < size_temp; i++)
    {
        if (PATHS[i].m_i[brenching_point] < 0)
        {
            for (int k = 0; k < PATHS_producing.size(); k++)
            {
                Merge_2_pathways(PATHS,
                                 PATHS[PATHS_producing[k]],
                                 PATHS[i],
                                 brenching_point,
                                 kinetics,
                                 rates);

                PATHS[PATHS_producing[k]].deleted = true;
            }
            PATHS[i].deleted = true;
        }
    }

    size_t temp = PATHS.size();
    for (size_t i = 0; i < temp; i++)
    {
        if (rates.delta_c(brenching_point) < 0 && PATHS[temp - 1 - i].m_i[brenching_point] < 0)
        {
            pathway P_n = PATHS[temp - 1 - i];
            P_n.deleted = false;
            P_n.f_k = P_n.f_k * abs(rates.p_i(brenching_point) - rates.d_i(brenching_point)) / rates.D_i(brenching_point);
            PATHS.push_back(P_n);
        }

        if (rates.delta_c(brenching_point) > 0 && PATHS[temp - 1 - i].m_i[brenching_point] > 0)
        {
            pathway P_n = PATHS[temp - 1 - i];
            P_n.deleted = false;
            P_n.f_k = P_n.f_k * abs(rates.p_i(brenching_point) - rates.d_i(brenching_point)) / rates.D_i(brenching_point);
            PATHS.push_back(P_n);
        }

        if (PATHS[temp - 1 - i].deleted)
        {
            if (PATHS[temp - 1 - i].m_i[brenching_point] > 0)
            {
                double fk = (PATHS[temp - 1 - i].f_k * dtilde
                             / rates.D_i(brenching_point));
                Update_deleted(PATHS[temp - 1 - i], rates, fk);
            }
            else
            {
                double fk = (PATHS[temp - 1 - i].f_k * ptilde
                             / rates.D_i(brenching_point));
                Update_deleted(PATHS[temp - 1 - i], rates, fk);
            }

            PATHS.erase(PATHS.begin() + temp - 1 - i);
        }
    }

    Update_deleted_for_branching_point(PATHS, rates, brenching_point);
}



// Merge 2 pathways
void Merge_2_pathways(vector<pathway> &PATHS,
                      const pathway   &Pathway_producing,
                      const pathway   &Pathway_consuming,
                      const int       &brenching_point,
                      const In_data   &kinetics,
                      Rates           &rates)
{
    pathway P_n;
    P_n.deleted = false;
    P_n.simple = false;

    // Merging paths from pathway_producing and pathway_consuming
    P_n.path.insert(P_n.path.end(), Pathway_producing.path.begin(),
                    Pathway_producing.path.end());
    P_n.path.insert(P_n.path.end(), Pathway_consuming.path.begin(),
                    Pathway_consuming.path.end());

    // Setting number of molecules m_i
    for (int i = 0; i < n_S; i++)
        P_n.m_i.push_back(Pathway_consuming.m_i[i] * Pathway_producing.m_i[brenching_point]
                          + Pathway_producing.m_i[i] * abs(Pathway_consuming.m_i[brenching_point]));

    doublearray1d x_jk(n_R);
    doublearray1d x_jl(n_R);
    doublearray1d x_jn(n_R);

    x_jk.setToValue(0);
    x_jl.setToValue(0);
    x_jn.setToValue(0);

    //creating vector x_jk
    for (size_t j = 0; j < Pathway_producing.path.size(); j++)
        x_jk(Pathway_producing.path[j]) = Pathway_producing.x_j[j];

    //creating vector x_jl
    for (size_t j = 0; j < Pathway_consuming.path.size(); j++)
        x_jl(Pathway_consuming.path[j]) = Pathway_consuming.x_j[j];

    //calculating x_jn with the formula (34)
    for (int j = 0; j < n_R; j++)
    {
        x_jn(j) = abs(Pathway_consuming.m_i[brenching_point]) * x_jk(j) +
                  Pathway_producing.m_i[brenching_point] * x_jl(j);
    }

    // adding non-zero terms to P_n
    for (size_t j = 0; j < Pathway_producing.path.size(); j++)
        P_n.x_j.push_back(x_jn(Pathway_producing.path[j]));

    for (size_t j = 0; j < Pathway_consuming.path.size(); j++)
        P_n.x_j.push_back(x_jn(Pathway_consuming.path[j]));

    // Calculating the P_n pathway rate
    P_n.f_k = ( Pathway_consuming.f_k * Pathway_producing.f_k ) / rates.D_i(brenching_point);
    Factorise(P_n);

    Sort_and_delete_pathway(P_n.path, P_n.x_j);

    PATHS.push_back(P_n);
}


// Updates the production (p) and destruction (d) arrays after merging
// pathways pw1 and pw2 into pwmerged.
void Calculate_rates(vector<pathway> &PATHS,
                     const In_data   &kinetics,
                     Rates           &rates)
{
    for (int i = 0; i < n_S; i++)
    {
        rates.p_i(i) = rates.tilda_p_i(i);
        rates.d_i(i) = rates.tilda_d_i(i);

        for (int k = 0; k < PATHS.size(); k++)
        {
            if (PATHS[k].m_i[i] > 0)
                rates.p_i(i) = rates.p_i(i) + PATHS[k].m_i[i] * PATHS[k].f_k;
            else
                rates.d_i(i) = rates.d_i(i) + abs(PATHS[k].m_i[i]) * PATHS[k].f_k;
        }

        rates.delta_i(i) = rates.p_i(i) - rates.d_i(i);

        rates.D_i(i) = max(rates.p_i(i), rates.d_i(i));
    }

}

// Elimination of insignificant pathways
void Insignificant(vector<pathway> &PATHS,
                   Rates           &rates,
                   const In_data   &input)
{
    Sort_PATH_by_rate(PATHS);

    for (int i = 0; i < PATHS.size(); i++)
    {
        if ((f_min > 0) && (PATHS[i].f_k < f_min))
        {
            PATHS[i].deleted = true;
            Update_deleted(PATHS[i], rates, PATHS[i].f_k);
        }
        else if ((max_path > 0) && (i > max_path))
        {
            PATHS[i].deleted = true;
            Update_deleted(PATHS[i], rates, PATHS[i].f_k);
        }
        //        } else if (PATHS[i].f_k == 0.0) {
        //            PATHS[i].deleted = true;
        //            Update_deleted(PATHS[i], rates, PATHS[i].f_k);
        //        }
    }

    size_t temp = PATHS.size();
    for (size_t i = 0; i < temp; i++)
        if (PATHS[temp - 1 - i].deleted) PATHS.erase(PATHS.begin() + temp - 1 - i);
    cout << "There are " << PATHS.size() << " pathways.\n" << endl;

}

// Sort by rate
void Sort_PATH_by_rate(vector<pathway> &PATHS)
{
    sort(PATHS.begin(), PATHS.end(), compare_pathway_rate);
}

bool compare_pathway_rate(const pathway &a,
                          const pathway &b)
{
    return (a.f_k > b.f_k);
}

// Updates the deleted rates after merging two pathways
void Update_deleted_after_merging(vector<pathway> &PATHS,
                                  pathway         &producing,
                                  pathway         &consuming,
                                  Rates           &rates,
                                  double           dtilde,
                                  double           ptilde,
                                  int              branching_point)
{
    double fk;

    // p[b] -> S[b] -> dtilde[b]
    fk = (producing.f_k
          * dtilde
          / rates.D_i(branching_point));
    Update_deleted(producing, rates, fk);


    // ptilde[b] -> S[b] -> d[b]
    fk = (consuming.f_k
          * ptilde
          / rates.D_i(branching_point));
    Update_deleted(consuming, rates, fk);
}


// Updates the contributions from deleted pathways, applying (43)-(45).
// This function is called both when we drop pathways in Insignificant
// and when we combine pathways through a branching point.
void Update_deleted(pathway         &p,
                    Rates           &rates,
                    double           f)
{
    for (int j = 0; j < p.path.size(); j++)
    {
        rates.tilda_r_j(p.path[j]) =
            rates.tilda_r_j(p.path[j]) + f * p.x_j[j];
    }

    for (int i = 0; i < n_S; i++)
    {
        if (p.m_i[i] > 0)
            rates.tilda_p_i(i) = rates.tilda_p_i(i) + p.m_i[i] * f;
        else
            rates.tilda_d_i(i) = rates.tilda_d_i(i) + abs(p.m_i[i]) * f;
    }
}

void Update_deleted_for_branching_point(vector<pathway> &PATHS,
                                        Rates           &rates,
                                        int              bp)
{
    double pd = rates.tilda_p_i(bp) * rates.tilda_d_i(bp) / rates.D_i(bp);

    rates.tilda_p_i(bp) -= pd;
    rates.tilda_d_i(bp) -= pd;
}


// Create elementary pathways
void Make_elementary(vector<pathway>    &PATHS,
                     vector<elementary> &Elementary,
                     const vector<int>  &all_brenching_points,
                     const Rates        &rates,
                     const In_data      &kinetics)
{
    for (int i = 0; i < PATHS.size(); i++)
    {
        if (!PATHS[i].simple)
        {
            Step_1(PATHS[i], Elementary[i], rates, kinetics);

            Add_pseudo_pathways(PATHS[i], Elementary[i], all_brenching_points);

            // Step 2 and Step 3
            Step_2(Elementary[i], all_brenching_points, kinetics);

            Remove_pseudo_pathways(Elementary[i]);
        }
    }
}


// Implements step_1 from 5.5.1
void Step_1(const pathway     &P_n,
            elementary        &elementary_paths,
            const Rates       &rates,
            const In_data     &kinetics)
{
    for (int j = 0; j < P_n.path.size(); j++)
    {
        pathway temp;
        temp.path.push_back(P_n.path[j]);
        temp.x_j.push_back(P_n.x_j[j]);
        for (int k = 0; k < n_S; k++)
        {
            temp.m_i.push_back(kinetics.matrix(k, P_n.path[j])*P_n.x_j[j]);
            temp.virtualm_i.push_back(0);
        }
        temp.deleted = false;

        temp.f_k = 0.0;

        elementary_paths.EL_pathways.push_back(temp);
    }
}

// Adds pseudo-pathways to elementaries for all branching points so far
void Add_pseudo_pathways(const pathway     &P_n,
                         elementary        &elementaries,
                         const vector<int> &bp)
{
    for (int i = 0; i < bp.size(); i++)
    {
        if (P_n.m_i[bp[i]] != 0)
        {
            pathway p;

            // HACK ALERT:
            // We will use reaction indices starting at n_R for the reactions in
            // the pseudo pathway.  This looks now like the best option since
            // A. Pathways are defined by their contained reactions
            // B. Setting an empty list makes (50) fail bc if the LHS is
            //    [], it will always evaluate to true.
            // C. If we add all reactions in P_n to the pseudo-pathway,
            //    (50) again gives problems bc it always evaluates to true
            //    whenever there is a pseudo-pw in the RHS.
            // D. Using negative indices is also a bad idea, since we will use
            //    them to access memory.
            // We must be careful that we do not do anything nasty with these
            // indices, and we must remember to remove them
            // below, in Remove_pseudo_pathways.
            p.path.push_back(n_R + bp[i]);
            p.x_j.push_back(abs(P_n.m_i[bp[i]]));

            p.m_i.assign(n_S, 0.0);
            p.virtualm_i.assign(n_S, 0.0);

            p.m_i[bp[i]] = -P_n.m_i[bp[i]];
            p.virtualm_i[bp[i]] = -P_n.m_i[bp[i]];

            // I am not sure if we need this
            p.f_k = 0.0;

            p.deleted = false;
            p.simple = true;
            elementaries.EL_pathways.push_back(p);
        }
    }
}


// Removes the pseudo-reactions from a pathway.
void Remove_pseudo_pathways(elementary &elementaries)
{
    size_t size = elementaries.EL_pathways.size();
    for (int i = 0; i < size; i++)
    {
        for (int k = 0; k < n_S; k++)
        {
            elementaries.EL_pathways[i].m_i[k] -= elementaries.EL_pathways[i].virtualm_i[k];
        }

        // The pseudo pathways contain negative reaction indices.  We have to
        // remove them.
        size_t kmax = elementaries.EL_pathways[i].path.size();

        for (int k = 0; k < kmax; k++)
        {
            // I love C++!! It's so concise...
            if (elementaries.EL_pathways[i].path[kmax - k - 1] >= n_R)
            {
                elementaries.EL_pathways[i].path.erase
                (elementaries.EL_pathways[i].path.begin() + kmax - k - 1);
                elementaries.EL_pathways[i].x_j.erase
                (elementaries.EL_pathways[i].x_j.begin() + kmax - k - 1);
            }
        }

        if (elementaries.EL_pathways[i].x_j.size() == 0)
        {
            printf ("Fuck off!");
            exit(EXIT_FAILURE);
        }
    }
}


// Implements step_2 from 5.5.1
void Step_2(elementary        &e_path,
            const vector<int> &all_brenching_points,
            const In_data     &kinetics)
{
    // loop over all brenching points befor index_bp
    for (int j = 0; j < all_brenching_points.size(); j++)
    {
        // Step 2.1,  P_tilda
        elementary P_tilda;

        // Step 2.2
        for (int k = 0; k < e_path.EL_pathways.size(); k ++)
        {
            if (e_path.EL_pathways[k].m_i[all_brenching_points[j]] == 0)
            {
                pathway tempo;
                Copy_pathway(e_path.EL_pathways[k], tempo);
                P_tilda.EL_pathways.push_back(tempo);
                Del_pathway(tempo);
            }
        }
        // Step 2.3
        Step_2_3(e_path, P_tilda, all_brenching_points[j], kinetics);

        // Step 2.4
        vector<pathway> ().swap(e_path.EL_pathways);
        for (int i = 0; i < P_tilda.EL_pathways.size(); i++)
            e_path.EL_pathways.push_back(P_tilda.EL_pathways[i]);

    }
}

// Implements step_2_3 from 5.5.1
void Step_2_3(const elementary &e_path,
              elementary       &P_tilda,
              int               bp,
              const In_data    &kinetics)
{
    vector<pathway> Production;
    vector<int> Production_index;

    for (int i = 0; i < e_path.EL_pathways.size(); i++)
    {
        if (e_path.EL_pathways[i].m_i[bp] > 0)
        {
            Production.push_back(e_path.EL_pathways[i]);
            Production_index.push_back(i);
        }
    }

    for (int i = 0; i < e_path.EL_pathways.size(); i++)
    {
        if (e_path.EL_pathways[i].m_i[bp] < 0)
        {
            for (int k = 0; k < Production.size(); k++)
            {
                pathway combination;

                Combine_2_pathways(combination, Production[k], e_path.EL_pathways[i], bp, kinetics, e_path);

                Condition_2_3(e_path, P_tilda, combination, Production[k],
                              e_path.EL_pathways[i],
                              Production_index[k], i,
                              kinetics);

            }

        }
    }

}

// Combine 2 pathways
void Combine_2_pathways(pathway           &combination,
                        const pathway     &Pathway_producing,
                        const pathway     &Pathway_consuming,
                        const int         &brenching_point,
                        const In_data     &input_raw,
                        const elementary  &e_path)
{
    combination.deleted = false;

    // Merging paths from pathway_producing and pathway_consuming
    combination.path.insert(combination.path.end(), Pathway_producing.path.begin(), Pathway_producing.path.end());
    combination.path.insert(combination.path.end(), Pathway_consuming.path.begin(), Pathway_consuming.path.end());

    // Setting number of molecules m_i
    for (int i = 0; i < n_S; i++)
    {
        combination.m_i.push_back(Pathway_consuming.m_i[i] * Pathway_producing.m_i[brenching_point]
                                  + Pathway_producing.m_i[i] * abs(Pathway_consuming.m_i[brenching_point]));

        combination.virtualm_i.push_back(Pathway_consuming.virtualm_i[i] * Pathway_producing.m_i[brenching_point]
                                         + Pathway_producing.virtualm_i[i] * abs(Pathway_consuming.m_i[brenching_point]));

    }
    doublearray1d x_jk(n_R + n_S);
    doublearray1d x_jl(n_R + n_S);
    doublearray1d x_jn(n_R + n_S);

    x_jk.setToValue(0);
    x_jl.setToValue(0);
    x_jn.setToValue(0);

    //creating vector x_jk
    for (size_t j = 0; j < Pathway_producing.path.size(); j++)
        x_jk(Pathway_producing.path[j]) = Pathway_producing.x_j[j];

    //creating vector x_jl
    for (size_t j = 0; j < Pathway_consuming.path.size(); j++)
    {
        x_jl(Pathway_consuming.path[j]) = Pathway_consuming.x_j[j];
    }

    for (int j = 0; j < n_R + n_S; j++)
    {
        x_jn(j) = abs(Pathway_consuming.m_i[brenching_point]) * x_jk(j) +
                  Pathway_producing.m_i[brenching_point] * x_jl(j);
    }

    // adding non-zero terms to P_n
    for (size_t j = 0; j < Pathway_producing.path.size(); j++)
        combination.x_j.push_back(x_jn(Pathway_producing.path[j]));

    for (size_t j = 0; j < Pathway_consuming.path.size(); j++)
        combination.x_j.push_back(x_jn(Pathway_consuming.path[j]));

    combination.f_k = 0.0;

    Sort_and_delete_pathway(combination.path, combination.x_j);

}


// Checks the validity of formula in 2.3
void Condition_2_3(const elementary &e_path,
                   elementary       &P_tilda,
                   const pathway    &combine,
                   const pathway     &Pathway_producing,
                   const pathway     &Pathway_consuming,
                   int k, int l,
                   const In_data    &kinetics)
{
    bool condition = true;
    for (int i = 0; i < e_path.EL_pathways.size(); i++)
    {
        if (e_path.EL_pathways[i].m_i.size() == 0)
        {
            cout << "Something is wrong in condition 2.3!" << endl;
            cout << e_path.EL_pathways.size() << endl;
            cout << i << endl;
            exit(1);
        }
        if ((i == k) || (i == l))
            continue;

        bool test = Contains(e_path.EL_pathways[i].path,
                             Pathway_consuming.path,
                             Pathway_producing.path);
        if (test)
        {
            condition = false;
            break;
        }

    }
    if (condition)
    {
        P_tilda.EL_pathways.push_back(combine);
    }
    else {  }
}


// Sort and delete repeted reactions in pathway paths
void Sort_and_delete_pathway(vector<int> &paths,
                             vector<int> &xj)
{
    int temp;

    // Sort elements in paths
    for (int i = 0; i < paths.size() - 1; i++)
    {
        for (int j = i + 1; j < paths.size(); j++)
            if (paths[i] > paths[j])
            {
                temp = paths[j];
                paths[j] = paths[i];
                paths[i] = temp;

                temp  = xj[j];
                xj[j] = xj[i];
                xj[i] = temp;
            }
    }

    // Delete repeted ones
    vector<int> del;
    for (int i = 0; i < paths.size() - 1; i++)
        if (paths[i] == paths[i + 1])
            del.push_back(i);

    for (int i = 0; i < del.size(); i++)
    {
        paths.erase(paths.begin() + del[del.size() - 1 - i]);
        xj.erase(xj.begin() + del[del.size() - 1 - i]);
    }

}


// Splitting into Sub-Pathways
void Splitting(vector<pathway>    &PATHS,
               vector<elementary> &Elementary,
               const vector<int>  &all_brenching_points,
               const Rates        &rates,
               const In_data      &kinetics)
{
    size_t max_size = PATHS.size();

    for (int i = 0; i < max_size; i++)
    {
        if (!PATHS[i].simple)
        {
            PATHS[i].deleted = false;
            doublearray1d w_k(Elementary[i].EL_pathways.size());

            for (int temp = 0; temp < Elementary[i].EL_pathways.size(); temp++)
                Factorise(Elementary[i].EL_pathways[temp]);

            if (w_k.getIndex1Size() == 1)
            {
                if (Paths_are_equal(Elementary[i].EL_pathways[0], PATHS[i]))
                {
                    PATHS[i].simple = true;
                }
                else
                {
                    cout << "Error: Problem with spliting." << endl;
                    Print_pathway(Elementary[i].EL_pathways[0], kinetics);
                    Print_pathway(PATHS[i], kinetics);
                    //exit (EXIT_FAILURE);
                }
            }
            else
            {

                Get_w_k(PATHS[i], Elementary[i], w_k);


                bool found = false;
                for (int k = 0; k < w_k.getIndex1Size(); k++)
                {
                    found = false;
                    for (int l = 0; l < max_size; l++)
                    {
                        if (Paths_are_equal(Elementary[i].EL_pathways[k], PATHS[l]))
                        {
                            PATHS[l].f_k = PATHS[l].f_k + w_k(k) * PATHS[i].f_k;
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        Elementary[i].EL_pathways[k].f_k = w_k(k) * PATHS[i].f_k;
                        Elementary[i].EL_pathways[k].deleted = false;
                        Elementary[i].EL_pathways[k].simple = true;
                        // add in PATH
                        Factorise( Elementary[i].EL_pathways[k]);
                        PATHS.push_back(Elementary[i].EL_pathways[k]);
                    }
                }
                PATHS[i].deleted = true;
            }
        }
    }

    // delete pathway from PATHS
    for (int i = 0; i < max_size; i++)
        if (PATHS[max_size - 1 - i].deleted)
            PATHS.erase(PATHS.begin() + max_size - 1 - i);
}


// Calculating w_k
void Get_w_k(pathway       &pw,
             elementary    &element,
             doublearray1d &w_k)
{
    double z = 0.0;
    doublearray1d x(pw.path.size());
    doublearray2d A(pw.path.size(), pw.path.size());
    doublearray1d b(pw.path.size());
    doublearray1d c(pw.path.size());

    Creat_LP(pw, element, A, b);
    for (int i = 0; i < pw.path.size(); i++) c(i) = 1;
    Solve_LP(z, x, A, b, c);
    for (int i = 0; i < w_k.getIndex1Size(); i++)
    {
        w_k(i) = x(i);
    }
}


// Creates square LP problem
void Creat_LP(const pathway    &pw,
              const elementary &element,
              doublearray2d    &A,
              doublearray1d    &b)
{
    doublearray2d A_0(n_R, element.EL_pathways.size());

    A_0.setToValue(0.0);
    for (int k = 0; k < element.EL_pathways.size(); k++)
        for (int j = 0; j < element.EL_pathways[k].x_j.size(); j++)
            A_0(element.EL_pathways[k].path[j], k) = element.EL_pathways[k].x_j[j];

    A.setToValue(0.0);

    /*
        Check if the problem is underdetermined.
    */
    // TODO: implement underdetermined case
    if (pw.path.size() < element.EL_pathways.size())
    {
        cout << endl;
        cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << endl;
        cout << "The LP problem is underdetermined. We suggest lowering the value of max_path." << endl;
        cout << endl;
        cout << "We will address this case in upcoming versions." << endl;
        cout << endl;
        cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << endl;
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < pw.path.size(); i++)
    {
        b(i) = pw.x_j[i];
        for (int k = 0; k < element.EL_pathways.size(); k++)
            A(i, k) = A_0(pw.path[i], k);
    }
}


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
              doublearray1d &c)
{
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_prob_name(lp, "Lp solver");
    glp_set_obj_dir(lp, GLP_MIN);

    // Create b
    glp_add_rows(lp, (int)b.getIndex1Size());
    for (int i = 0; i < b.getIndex1Size(); i++)
    {
        glp_set_row_name(lp, i + 1, "b");
        glp_set_row_bnds(lp, i + 1, GLP_FX, b(i), b(i));
    }

    // Create x
    glp_add_cols(lp, (int)x.getIndex1Size());
    for (int i = 0; i < x.getIndex1Size(); i++)
    {
        glp_set_col_name(lp, i + 1, "x");
        glp_set_col_bnds(lp, i + 1, GLP_LO, 0.0, 0.0); // x_i >= 0
        glp_set_obj_coef(lp, i + 1, c(i));
    }

    // The number of elements in the matrix A
    int size = (int)A.getIndex1Size();
    int ne = (int)pow(size, 2.0);

    int *ia = new int[ne + 1];
    int *ja = new int[ne + 1];
    double *ar = new double[ne + 1];

    for (int i = 1; i < size + 1; i++)
    {
        ia[i] = 1; ja[i] = i; ar[i] = A(0, i - 1);
    }


    int temp = size + 1;
    for (int j = 0; j < size; j++)
    {

        for (int i = 1; i < size; i++)
        {
            ia[temp] = i + 1;
            ja[temp] = j + 1;
            ar[temp] = A(i, j);

            temp = temp + 1;
        }
    }

    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;

    glp_load_matrix(lp, ne, ia, ja, ar);
    if (glp_simplex(lp, &parm))
    {
        cout << "LP problem can not be solved!" << endl;
        exit (EXIT_FAILURE);
    }

    z = glp_get_obj_val(lp);
    for (int i = 0; i < x.getIndex1Size(); i++)
        x(i) = glp_get_col_prim(lp, i + 1);

    glp_delete_prob(lp);

    delete []ia;
    delete []ja;
    delete []ar;
}


// GCD of x_j with f_k update
void Factorise (pathway &pw)
{
    int gcd = GCD(pw.x_j);
    if (gcd != 1)
    {
        pw.f_k = pw.f_k * gcd;
        for (int j = 0; j < pw.x_j.size(); j++)
            pw.x_j[j] = pw.x_j[j] / gcd;

        for (int i = 0; i < pw.m_i.size(); i++)
            pw.m_i[i] = pw.m_i[i] / gcd;

    }
}


// Copy pathway
void Copy_pathway(const pathway &from,
                  pathway       &to)
{
    Del_pathway(to);
    to.deleted = from.deleted;
    to.f_k = from.f_k;

    for (int i = 0; i < from.path.size(); i++)
    {
        to.path.push_back(from.path[i]);
        to.x_j.push_back(from.x_j[i]);
    }
    for (int i = 0; i < from.m_i.size(); i++)
    {
        to.m_i.push_back(from.m_i[i]);
        to.virtualm_i.push_back(from.virtualm_i[i]);
    }
}


// delete pathway
void Del_pathway(pathway &del)
{
    del.f_k = 0.0;
    del.deleted = false;
    vector<int> ().swap(del.m_i);
    vector<int> ().swap(del.virtualm_i);
    vector<int> ().swap(del.path);
    vector<int> ().swap(del.x_j);
}


// Checks if one \in in_two + in_three
bool Contains(const vector<int> &one,
              const vector<int> &in_two,
              const vector<int> &in_three)
{
    vector<int> second;
    second.insert(second.end(), in_two.begin(), in_two.end());
    second.insert(second.end(), in_three.begin(), in_three.end());

    for (int i = 0; i < one.size(); i++)
    {
        bool is_contained = false;
        for (int j = 0; j < second.size(); j++)
            if (one[i] == second[j])
            {
                is_contained = true;
                break;
            }

        if (!is_contained)
        {
            /* There is one element in one not contained in two + three. */
            return false;
        }
    }
    return true;
}


// Checking if pathway_1 = pathway_2
bool Paths_are_equal(pathway pathway_1,
                     pathway pathway_2)
{
    // TODO: Check also that x_j are equal
    if (pathway_1.path.size() != pathway_2.path.size()) return false;

    for (int i = 0; i < pathway_1.path.size() ; i++)
    {
        bool inside = false;
        for (int j = 0; j < pathway_2.path.size(); j++)
            if (pathway_1.path[i] == pathway_2.path[j])
                if (pathway_1.x_j[i] == pathway_2.x_j[j])
                    inside = true;
        if (inside == false) return false;
    }
    return true;
}

// Returns GCD
int GCD (vector<int> x_j)
{
    return std::accumulate(x_j.begin() + 1, x_j.end(), x_j[0], gcd_local);
    //accumulate(x_j.begin(), x_j.size(), 0);
}


int gcd_local(int a,
              int b)
{
    for (;;)
    {
        if (a == 0) return b;
        b %= a;
        if (b == 0) return a;
        a %= b;
    }
}












