#include "DispComb.h"
#include <cstddef> 

unsigned int* generateCombinations(unsigned int* arr, unsigned int n, unsigned int k, uint64_t numCombinations) {
    // Allocate memory for the results (numCombinations * k elements in total)
    unsigned int* result = (unsigned int*)malloc(numCombinations * k * sizeof(unsigned int));
    if (result == nullptr) {
        return nullptr;
    }

    // Array to store current combination
    unsigned int* indices = (unsigned int*)malloc(k * sizeof(unsigned int));
    if (indices == nullptr) {
        free(result);
        return nullptr;
    }

    // Initialize the first combination (the first k elements)
    for (unsigned int i = 0; i < k; ++i) {
        indices[i] = i;
    }

    uint64_t count = 0;

    // Generate all combinations
    while (count < numCombinations) {
        // Store the current combination into the result array
        for (unsigned int i = 0; i < k; ++i) {
            result[count * k + i] = arr[indices[i]];
        }
        count++;

        // Find the rightmost index that can be incremented
        int i = k - 1;
        while (i >= 0 && indices[i] == n - k + i) {
            --i;
        }

        // If no valid index is found, all combinations are generated
        if (i < 0) break;

        // Increment this index
        ++indices[i];

        // Update subsequent indices to maintain the lexicographical order
        for (unsigned int j = i + 1; j < k; ++j) {
            indices[j] = indices[j - 1] + 1;
        }
    }

    // Clean up dynamically allocated memory
    // free(indices);

    // Return the result array
    return result;
}


unsigned int* generateDispositions(unsigned int* arr, unsigned int n, unsigned int k, uint64_t numResults) {
    // Allocate memory for the results (numResults * k elements in total)
    unsigned int* result = (unsigned int*)malloc(numResults * k * sizeof(unsigned int));
    if (result == nullptr) {
        return nullptr;
    }

    // Array to keep track of the indices used in current disposition
    unsigned int* indices = (unsigned int*)malloc(k * sizeof(unsigned int));
    if (indices == nullptr) {
        // free(result);
        return nullptr;
    }

    // Array to track which elements are used
    bool* used = (bool*)malloc(n * sizeof(bool));
    if (used == nullptr) {
        free(result);
        free(indices);
        return nullptr;
    }

    // Initialize the used array to track available elements
    for (unsigned int i = 0; i < n; ++i) {
        used[i] = false;
    }

    // Initialize the first disposition to the first k elements
    for (unsigned int i = 0; i < k; ++i) {
        indices[i] = i;
        used[i] = true;
    }

    unsigned int count = 0;
    bool done = false;  // Flag to control the termination of the outer loop

    // Generate all possible dispositions
    while (!done && count < numResults) {
        // Copy the current disposition into the result array
        for (unsigned int i = 0; i < k; ++i) {
            result[count * k + i] = arr[indices[i]];
        }
        count++;

        // Start from the last element of indices and generate the next permutation
        int i = k - 1;
        bool foundNext = false;  // Flag to track if we found the next valid element

        // Traverse backward to find the next valid element
        while (i >= 0 && !foundNext) {
            used[indices[i]] = false;  // Mark current element as unused
            indices[i]++;  // Move to the next element

            // Try to find the next unused element
            while (indices[i] < n && used[indices[i]]) {
                indices[i]++;
            }

            if (indices[i] < n) {
                used[indices[i]] = true;  // Mark as used
                foundNext = true;  // We found the next valid element
            } else {
                i--;  // If no valid element is found, move back to the previous index
            }
        }

        // If no valid next element was found and i < 0, we are done
        if (!foundNext && i < 0) {
            done = true;  // Set done flag to terminate the loop
        } else if (foundNext) {
            // Reset the following indices and mark them as unused
            for (unsigned int j = i + 1; j < k; ++j) {
                indices[j] = 0;
                while (used[indices[j]]) {
                    indices[j]++;
                }
                used[indices[j]] = true;
            }
        }
    }

    // Clean up dynamically allocated memory
    free(used);
    free(indices);

    return result;
}

uint64_t binomialCoefficient(unsigned int n, unsigned int k) 
{
    if (k > n - k) // Take advantage of symmetry
        k = n - k;

    uint64_t result = 1;
    for (unsigned int i = 0; i < k; ++i) 
    {
        result *= (n - i);
        result /= (i + 1);
    }

    return result;
}

uint64_t permutation(unsigned int n, unsigned int k) 
{
    uint64_t result = 1;
    for (unsigned int i = 0; i < k; ++i) 
    {
        result *= n - i;
    }
    return result;
}


