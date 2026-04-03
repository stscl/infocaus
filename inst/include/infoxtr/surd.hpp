

#ifndef INFOXTR_SURD_HPP
#define INFOXTR_SURD_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <thread>
#include <cstdint>
#include <algorithm>
#include <stdexcept>

namespace infoxtr 
{

namespace surd 
{

    using DiscMat = std::vector<std::vector<uint64_t>>;
    using ContMat = std::vector<std::vector<double>>;

    /***********************************************************
     * Result structure
     ***********************************************************/
    struct SURDRes
    {
        std::vector<double> values;
        std::vector<uint8_t> types;
        std::vector<std::vector<size_t>> var_indices;

        size_t size() const noexcept { return values.size(); }
    };

    /***************************************************************
     * Synergistic-Unique-Redundant Decomposition for Discrete Data
     ***************************************************************/
    inline SURDRes surd(
        const DiscMat& mat,
        size_t max_order = std::numeric_limits<size_t>::max(),
        size_t threads = 1,
        double base = 2.0,
        bool normalize = false)

    /*****************************************************************
     * Synergistic-Unique-Redundant Decomposition for Continuous Data
     *****************************************************************/
    inline SURDRes surd(
        const ContMat& mat,
        size_t max_order = std::numeric_limits<size_t>::max(),
        size_t k = 3,
        size_t alg = 0,
        size_t threads = 1,
        double base = 2.0,
        bool normalize = false)
    
} // namespace surd

}

#endif // INFOXTR_SURD_HPP