/*
code developed and maintained by (jmw@ruc.edu.cn, RUC, China) date 2024
*/

#pragma once

#include <unordered_map>
#include <vector>
#include <string>

#include <bitset>
#include <stdexcept>

class AsnciSpace {
private:
    // Mapping from index to slater determinants (SDs), where the SDs contains two Ints
    std::unordered_map<int, std::vector<int>> slaterdets;
    int current_index = 0; // Current index, starting from 0

    std::vector<std::vector<int>> import_SDs;
public:
    // Constructor to initialize AsnciSpace with import_SDs
    AsnciSpace(const std::vector<std::vector<int>>& import_SDs);

    // Add an SD and return the corresponding index
    int addSD(const std::vector<int>& arr);

    // Get the SD by index
    std::vector<int> getSDByIndex(int index);

    // Get the index by SD
    int getIndexBySD(const std::vector<int>& SD);

    // Generate a new SD from important SDs
    std::vector<int> generate_new_SD(std::vector<std::vector<int>> import_SDs);

    // Derive the corresponding binary string from the SD
    std::string getBinaryStringFromSD(const std::vector<int>& SD);
};