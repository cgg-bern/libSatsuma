//  SPDX-FileCopyrightText: 2023 Martin Heistermann <martin.heistermann@unibe.ch>
//  SPDX-License-Identifier: MIT
#include <libsatsuma/Problems/BiMDF.hh>
#include <libsatsuma/Extra/Highlevel.hh>
#include <libsatsuma/IO/read_bimdf.hh>
#include <libTimekeeper/StopWatchPrinting.hh>
#include <map>


int main(int argc, char *argv[])
{
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <input.bimdf>\n" << std::endl;
        return 1;
    }
    using BiMDF = Satsuma::BiMDF;
    auto bimdf = Satsuma::read_bimdf(argv[1]);

    auto config = Satsuma::BiMDFSolverConfig {
        .matching_solver = Satsuma::MatchingSolver::Lemon
    };
    auto result = Satsuma::solve_bimdf(*bimdf, config);
    std::cout << "Total cost: " << result.cost << std::endl;

    std::cout << result.stopwatch << std::endl;
    return 0;
}
