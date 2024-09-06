#ifndef TIMER_CPP
#define TIMER_CPP
#include <iostream>
#include <chrono>
#include <string>
#include <sstream>
#include <iomanip>

class Timer {
public:
    Timer() : start_time_point(std::chrono::high_resolution_clock::now()) {}

    void reset() {
        start_time_point = std::chrono::high_resolution_clock::now();
    }

    double elapsedMilliseconds() const {
        auto end_time_point = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(start_time_point).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(end_time_point).time_since_epoch().count();
        return static_cast<double>(end - start) * 0.001; // Convert microseconds to milliseconds
    }

    double elapsedSeconds() const {
        return elapsedMilliseconds() / 1000.0; // Convert milliseconds to seconds
    }

    void printElapsed(const std::string& message = "", bool log = false) 
    {
        auto milliseconds = elapsedMilliseconds();
        auto seconds = elapsedSeconds();
        reset();

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2);
        oss << message << ". Elapsed time: " << milliseconds << " ms (" << seconds << " seconds)";

        if (log) 
            std::cout << oss.str() << "\n";
        //return oss.str();
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_point;
};


/*
Timer timer; // The timer starts automatically

// code to measure goes here

timer.printElapsed("First operation"); // Print elapsed time with a message

timer.reset(); // Optionally reset the timer for another measurement

// More code to measure goes here

timer.printElapsed("Second operation"); // Print elapsed time with a different message
*/

#endif