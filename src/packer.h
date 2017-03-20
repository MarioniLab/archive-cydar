#ifndef PACKER_H
#define PACKER_H
#include "cydar.h"

template<class Iter>
void pack_index_vector(std::deque<int>& sorted_ids, Iter first, Iter last) {
    std::sort(first, last);
    sorted_ids.clear();
    if (first==last) { return; }
    if (*first < 0) {
        throw std::runtime_error("all indices must be non-negative integers");
    }

    sorted_ids.push_back((*first)+1);
    Iter current=first, previous=first;

    while (1) { 
        ++current;

        // Looping across consecutive indices.
        // Also scans across repeated indices; array is sorted, so current<previous+1 only occurs with repeats.
        while (current!=last && (*current)<=(*previous)+1) {
            ++current;
            ++previous;
        }

        // If it's not the same as the last entry we added, the loop above must have done something.
        if ((*previous)+1!=sorted_ids.back()) {
            sorted_ids.push_back(-(*previous)-1);
        }   

        // Adding the next non-consecutive point.
        if (current==last) {
            break;
        } else {
            sorted_ids.push_back((*current)+1);
        }

        ++previous;
    }
    return;
}

template<class Iter>
void unpack_index_vector (std::deque<int>& output, Iter start, Iter end) {
    output.clear();
    while (start!=end) {
        if ((*start) > 0) {
            if (!output.empty() && (*start) < output.back()) {
                throw std::runtime_error("absolute values of compressed indices must always increase");
            }
            output.push_back(*start);
        } else if ((*start) < 0) {
            if (output.empty() || *(start-1) < 0 || output.back() > -(*start)) {
                throw std::runtime_error("inappropriate negative values in compressed index vector");
            }
            while (output.back()!=-(*start)) {
                output.push_back(output.back()+1);
            }
        } else {
            throw std::runtime_error("zero values in compressed index vector");
        }
        ++start;
    }
    return;
}

#endif
