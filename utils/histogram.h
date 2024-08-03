

// S
#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <string>
#include <boost/format.hpp>

#include "error.h"


template <typename T>
class Histogram {

    // Define upper and lower bounds
    const T lower;
    const T upper;
    const int bins;



    // Bucket struct
    struct Bucket {
        std::string mark;
        T count;
        double frequency;
    };

    // Containers to hold the buckets and counts
    // The number of buckets is bins + 2 because there
    // is one extra bucket on each side to hold values
    // that are outside the range
    std::vector<Bucket> buckets;

    // Bucket size
    const T bs;

    // total number of elements inserted
    int samples;
    T max;
    T min;
    T sum;


  public:

    Histogram(T lower, T upper, int bins)
        : lower(lower), upper(upper), bins(bins),
          bs((upper - lower) / static_cast<T>(bins)),
          samples(0), max(0), min(0), sum(0)
    {

        panic_if(bs <= 0, "Bucket size must be positive");

        buckets.resize(bins + 2);
        buckets[0] = Bucket{
            boost::str(boost::format("< %.1f") % lower), 0, 0};
        buckets[bins + 1] = Bucket{
            boost::str(boost::format("> %.1f") % upper), 0, 0};

        for (int i = 1; i <= bins; ++i) {
            if (bs == 1) {
                buckets[i] = Bucket{
                    boost::str(boost::format("%d") %
                        (lower + bs * static_cast<T>(i-1))),
                    0,
                    0
                };
            } else {
                buckets[i] = Bucket{
                    boost::str(boost::format("%.1f-%.1f") %
                        (lower + bs * static_cast<T>(i-1)) %
                        (lower + bs * static_cast<T>(i))),
                    0,
                    0
                };
            }
        }

    }

    void insert(T value, int count=1) {
        samples += count;
        sum += value * count;
        if (value > max) {
            max = value;
        }
        if (value < min) {
            min = value;
        }

        if (value < lower) {
            buckets[0].count += count;
            return;
        } else if (value > upper) {
            buckets[bins + 1].count += count;
            return;
        }

        int idx = static_cast<int>((value - lower) / bs) + 1;
        buckets[idx].count += count;
    }

    void insert(const std::vector<T>& values) {
        for (auto v : values) {
            insert(v);
        }
    }

    void reset() {
        samples = 0;
        max = 0;
        min = 0;
        sum = 0;
        for (auto& b : buckets) {
            b.count = 0;
        }
    }


    std::string print(bool perc=false, bool cdf=false, int width = 40) {
        T _max = 0;
        for (const auto& b : buckets) {
            if (b.count > _max) {
                _max = b.count;
            }
        }

        std::ostringstream res;
        res << boost::format("N:%d Min:%d, Max:%d, Sum:%d, Avg:%.2f\n")
            % samples % min % max % sum
            % getAvg();


        int cum = 0;
        for (const auto& b : buckets) {
            int barLen = _max > 0 ? static_cast<int>((b.count * width + _max / 2) / _max) : 0;

            res << boost::format("%-10s [%-4d]\t") % b.mark % b.count;

            if (perc) {
                res << boost::format("%.1f%%\t") % (100.0 * b.count / samples);
            }

            if (cdf) {
                cum += b.count;
                res << boost::format("%.1f%%\t") % (100.0 * cum / samples);
            }

            res << boost::format("|%s\n") % std::string(barLen, '*');

        }
        return res.str();
    }

    std::string printCDF() {
        return print(false, true);
    }

    double getAvg() {
        return static_cast<double>(sum) / static_cast<double>(samples);
    }

    T getMax() {
        return max;
    }

    T getMin() {
        return min;
    }

};

