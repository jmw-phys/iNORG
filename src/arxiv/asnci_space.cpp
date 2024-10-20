#include "asnci_space.h"
// Constructor to initialize AsnciSpace with import_SDs
AsnciSpace::AsnciSpace(const std::vector<std::vector<int>>& import_SDs) {
    // Add the initial import_SDs to slaterdets
    for (const auto& sd : import_SDs) {
        addSD(sd);
    }

    // Generate new SDs by applying H and add them to slaterdets
    for (const auto& sd : import_SDs) {
        std::vector<int> new_sd = generate_new_SD({sd});
        addSD(new_sd);
    }
}

int AsnciSpace::addSD(const std::vector<int>& arr) {
    if (arr.size() != 2) {
        throw std::invalid_argument("Array must contain exactly 2 elements.");
    }
    slaterdets[current_index] = arr;
    return current_index++;
}

std::vector<int> AsnciSpace::getSDByIndex(int index) {
    if (slaterdets.find(index) != slaterdets.end()) {
        return slaterdets[index];
    } else {
        return {}; // Return an empty vector if not found
    }
}

int AsnciSpace::getIndexBySD(const std::vector<int>& SD) {
    for (const auto& pair : slaterdets) {
        if (pair.second == SD) {     // ! Here need to change!
            return pair.first;
        }
    }
    return -1; // Return -1 if not found
}

std::vector<int> AsnciSpace::generate_new_SD(std::vector<std::vector<int>> import_SDs) {
    // ! Here need to change!
    // Example logic to generate a new SD from import_SDs
    // Here we simply add the first elements of all import_SDs
    int new_sd_1 = 0;
    int new_sd_2 = 0;
    for (const auto& sd : import_SDs) {
        if (sd.size() == 2) {
            new_sd_1 += sd[0];
            new_sd_2 += sd[1];
        }
    }
    return {new_sd_1, new_sd_2};
}

std::string AsnciSpace::getBinaryStringFromSD(const std::vector<int>& SD) {
    if (SD.size() != 2) {
        throw std::invalid_argument("Array must contain exactly 2 elements.");
    }
    std::bitset<32> bits1(SD[0]);
    std::bitset<32> bits2(SD[1]);
    return bits1.to_string() + bits2.to_string();
}