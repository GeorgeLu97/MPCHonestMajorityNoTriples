#include "Timer.h"

#include <stdlib.h>
#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>
#include <chrono>
#include <libscapi/include/infra/Measurement.hpp>

using namespace std;
using namespace std::chrono;


Timer::Timer(string & outputFile,
             unsigned int numberOfIterations)
  : _outputFile(outputFile), _numberOfIterations(numberOfIterations) {}

int Timer::getTaskIdx(string name)
{
    auto it = find(_names.begin(), _names.end(), name);
    auto idx = distance(_names.begin(), it);
    if (it == _names.end()) {
        _names.push_back(name);
        _startTimes.push_back(vector<double>(_numberOfIterations, 0));
        _totalTimes.push_back(vector<double>(_numberOfIterations, 0));
        _isTiming.push_back(vector<bool>(_numberOfIterations, false));
    }
    return idx;
}

void Timer::startSubTask(string taskName, unsigned int currentIterationNum)
{
    if (currentIterationNum >= _numberOfIterations) {
        throw runtime_error("Current iteration "
                            + to_string(currentIterationNum)
                            + " is out of range");
    }

    int taskIdx = getTaskIdx(taskName);
    if (_isTiming[taskIdx][currentIterationNum]) {
        throw runtime_error("Already start a timing for " + taskName
                            + "(" + to_string(currentIterationNum) + ")");
    }

    auto now = system_clock::now();
    auto ms = (double) time_point_cast<nanoseconds>(now).time_since_epoch()
                                                        .count() / 1000000;
    _startTimes[taskIdx][currentIterationNum] = ms;
    _isTiming[taskIdx][currentIterationNum] = true;
}

void Timer::endSubTask(string taskName, unsigned int currentIterationNum)
{
    if (currentIterationNum >= _numberOfIterations) {
        throw runtime_error("Current iteration "
                            + to_string(currentIterationNum)
                            + " is out of range");
    }

    int taskIdx = getTaskIdx(taskName);
    if (!_isTiming[taskIdx][currentIterationNum]) {
        throw runtime_error("No timing for " + taskName
                            + "(" + to_string(currentIterationNum) + ")");
    }

    auto now = system_clock::now();
    auto ms = (double) time_point_cast<nanoseconds>(now).time_since_epoch()
                                                        .count() / 1000000;
    _totalTimes[taskIdx][currentIterationNum] += ms
                      - _startTimes[taskIdx][currentIterationNum];
    _isTiming[taskIdx][currentIterationNum] = false;
}

void Timer::analyze()
{
    char buff[255];
    string filePath = getcwd(buff, 255);
    string fileName = filePath + "/" + _outputFile + ".json";

    //party is the root of the json objects
    json party = json::array();

    for (int taskNameIdx = 0; taskNameIdx < _names.size(); taskNameIdx++)
    {
        //Write for each task name all the iteration
        json task = json::object();
        task["name"] = _names[taskNameIdx];

        for (int iterationIdx = 0; iterationIdx < _numberOfIterations; iterationIdx++)
        {
            ostringstream streamObj;
            streamObj << fixed << setprecision(3)
                      << _totalTimes[taskNameIdx][iterationIdx];
            task["iteration_" + to_string(iterationIdx)] = streamObj.str();
        }

        party.insert(party.begin(), task);
    }

    //send json object to create file
    try {
        ofstream myfile (fileName, ostream::out);
        myfile << party;
    } catch (exception& e) {
        cout << "Exception thrown : " << e.what() << endl;
    }
}

Timer::~Timer()
{
    analyze();
}

