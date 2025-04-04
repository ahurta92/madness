//
// Created by Florian Bischoff on 2/15/24.
//

#ifndef MP3_H
#define MP3_H


#include <madness/mra/mra.h>
#include<madness/mra/commandlineparser.h>
#include<madness/chem/ccpairfunction.h>
#include<madness/chem/CCStructures.h>
#include<madness/chem/CCPotentials.h>
#include<madness/mra/QCCalculationParametersBase.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <madness/mra/macrotaskq.h>

namespace madness {
class MP3 {
public:

    // MP3(World& world, const std::shared_ptr<Nemo> nemo, const CCParameters& param) {}

    MP3(World& world, const Info& info) : world(world), info(info) {}

    double mp3_energy_contribution(const Pairs<CCPair>& mp2pairs) const;

    /// compute the MP3 energy contribution, macrotask version
    double mp3_energy_contribution_macrotask_driver(const Pairs<CCPair>& mp2pairs) const;

private:
    World& world;
    const Info& info;

    /// helper class for calculating the MP3 energy contributions
    class MacroTaskMP3 : public MacroTaskOperationBase {

        class Partitioner : public MacroTaskPartitioner {
        public:
            Partitioner(const std::string shape) {
                min_batch_size=1;
                max_batch_size=1;
                if (shape=="triangular") dimension=1;
                else if (shape=="square") dimension=2;
                else {
                    std::string msg = "Unknown partitioning shape: " + shape;
                    MADNESS_EXCEPTION(msg.c_str(), 1);
                }
            };
        };

    public:
        MacroTaskMP3(const std::string diagram) : diagram(diagram) {
            name="MP3_"+diagram;
            std::string shape;
            if (diagram=="ghij" or diagram=="cd" or diagram=="ef") shape="triangular";
            else if (diagram=="klmn") shape="square";
            else {
                std::string msg = "Unknown MP3 diagram: " + diagram;
                MADNESS_EXCEPTION(msg.c_str(), 1);
            }
            partitioner.reset(new Partitioner(shape));
        }
        std::string diagram="unknown";

        typedef std::tuple<
                const std::vector<int>&,        // dummy vector of size npair or nocc for scheduling
                const std::vector<int>&,        // dummy vector of size npair or nocc for scheduling
                const std::vector<std::vector<CCPairFunction<double,6>>>& ,                 // all pairs ij
                const Info&,
                const std::vector<std::string>& > argtupleT;

        using resultT =ScalarResult<double>;

        resultT allocator(World& world, const argtupleT& argtuple) const {
            return ScalarResult<double>(world);
        }

        resultT operator() (const std::vector<int>& ij_vec,                         // dummy vector of size npair or nocc
                            const std::vector<int>& j_vec,                          // dummy vector of size 0 or nocc
                            const std::vector<std::vector<CCPairFunction<double,6>>>& pair_square,                 // all pairs ij
                            const Info& info,
                            const std::vector<std::string>& argument) const {

            // the partitioner will break the input vector of pairs into single pairs
            MADNESS_CHECK(ij_vec.size()==1);
            MADNESS_CHECK(batch.result.size()==1);

            // nact=active occupied orbitals
            const long nact=info.mo_ket.size()-info.parameters.freeze();
            MADNESS_CHECK(pair_square.size()==nact*nact);

            // loop over pairs i<j (triangular) or i,j (square)
            bool is_triangular=(j_vec.size()==0);
            if (is_triangular) MADNESS_CHECK(partitioner->dimension==1);

            // determine the orbital indices i and j for the pair
            int i=0, j=0;
            if (is_triangular) {
                // the batch index is the ij composite index [0,nact*(nact+1)-1]
                const long ij=batch.result.begin;
                // turn composite index ij into i and j, taking care of frozen orbitals
                PairVectorMap tri_map=PairVectorMap::triangular_map(info.parameters.freeze(),info.mo_ket.size());
                auto ij_to_i_and_j = [&tri_map](const int ij) { return tri_map.map[ij]; };
                auto [ii,jj]=ij_to_i_and_j(ij);
                i=ii;
                j=jj;
            } else {
                MADNESS_CHECK(partitioner->dimension==2);
                MADNESS_CHECK(j_vec.size()==ij_vec.size());
                i=batch.input[0].begin+info.parameters.freeze();
                j=batch.input[1].begin+info.parameters.freeze();
            }
            // print("i,j,parameters.freeze()=",i,j,parameters.freeze());

            // convert vector of vectors back to Pairs
            PairVectorMap square_map=PairVectorMap::quadratic_map(info.parameters.freeze(),info.mo_ket.size());
            auto clusterfunctions=Pairs<std::vector<CCPairFunction<double,6>>>::vector2pairs(pair_square,square_map);

            double result=0.0;
            World& world=info.R_square.world();
            if (diagram=="cd")
                result= MP3::compute_mp3_cd(world,i,j,clusterfunctions,info,argument);
            else if (diagram=="ef")
                result= MP3::compute_mp3_ef(world,i,j,clusterfunctions,info,argument);
            else if (diagram=="ghij")
                result= MP3::compute_mp3_ghij(world,i,j,clusterfunctions,info,argument);
            else if (diagram=="klmn")
                result= MP3::compute_mp3_klmn(world,i,j,clusterfunctions,info,argument);
            else {
                std::string msg = "Unknown MP3 diagram: " + diagram;
                MADNESS_EXCEPTION(msg.c_str(), 1);
            }
            auto result1=ScalarResult<double>(world);
            result1=result;
            return result1;

        };


    };


    double compute_mp3_cd(const Pairs<CCPair>& mp2pairs, const Info& info) const;
    double compute_mp3_ef(const Pairs<CCPair>& mp2pairs, const Info& info) const;
    double compute_mp3_ef_with_permutational_symmetry(const Pairs<CCPair>& mp2pairs, const Info& info) const;
    double compute_mp3_ef_low_scaling(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions, const Info& info) const;
    double compute_mp3_ef_as_overlap(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions, const Info& info) const;
    double compute_mp3_ghij(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions, const Info& info) const;
    double compute_mp3_ghij_fast(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions, const Info& info) const;
    double compute_mp3_klmn(const Pairs<CCPair>& mp2pairs, const Info& info) const;
    double compute_mp3_klmn_fast(const Pairs<CCPair>& mp2pairs, const Info& info) const;
    double mp3_test(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions, const Info& info) const;

    /// compute the cd term for single pair ij
    static double compute_mp3_cd(World& world,
                                 const long i, const long j,
                                 const Pairs<std::vector<CCPairFunction<double,6>>>& pair_square,
                                 const Info& info,
                                 const std::vector<std::string>& argument);

    /// compute the ef term for single pair ij
    static double compute_mp3_ef(World& world,
                                 const long i, const long j,
                                 const Pairs<std::vector<CCPairFunction<double,6>>>& pair_square,
                                 const Info& info,
                                 const std::vector<std::string>& argument);

    /// compute the ghij term for single pair ij
    ///
    /// the term actually scales linearly with the number of occupied orbitals i, so for all i!=j return zero
    static double compute_mp3_ghij(World& world,
                                   const long i, const long j,
                                   const Pairs<std::vector<CCPairFunction<double,6>>>& pair_square,
                                   const Info& info,
                                   const std::vector<std::string>& argument);

    /// compute the klmn term for single pair ij
    static double compute_mp3_klmn(World& world,
                                   const long i, const long j,
                                   const Pairs<std::vector<CCPairFunction<double,6>>>& pair_square,
                                   const Info& info,
                                   const std::vector<std::string>& argument);

};
}


#endif //MP3_H
