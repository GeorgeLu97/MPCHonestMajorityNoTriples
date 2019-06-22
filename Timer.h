#ifndef SP_TIMER_H_
#define SP_TIMER_H_

#include <stdlib.h>
#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>
#include <chrono>
#include <libscapi/include/infra/Measurement.hpp>

using namespace std;
using namespace std::chrono;

class Timer
{
private:
    string _outputFile;
    unsigned int _numberOfIterations;
    vector<string> _names;
    vector<vector<bool>> _isTiming;
    vector<vector<double>> _startTimes;
    vector<vector<double>> _totalTimes;

public:
    Timer(string & _outputFile, unsigned int numberOfIterations);
    ~Timer();

    void startSubTask(string taskName, unsigned int currentIterationNum);
    void endSubTask(string taskName, unsigned int currentIterationNum);
    void analyze();

private:
    int getTaskIdx(string name);
};

#endif /* TIMER_H_ */
