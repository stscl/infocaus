

#ifndef INFOXTR_SURD_HPP
#define INFOXTR_SURD_HPP

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

     /****************************************************************
     * Synergistic-Unique-Redundant Decomposition for Continuous Data
     *****************************************************************/
    
} // namespace surd

}


#endif // INFOXTR_SURD_HPP