#include <iostream>
#include "MyStatistics.h"

using std::cerr;

int main(){
    MyStatistics<int, 5> data;
    MyStatistics<int, 5>::newArray na = {2, 5, 4, 23, 12};
    data.addRecord(na);
    data.addRecord({16, 44, 20, 51, 22});
    data.print();
    cerr << "record count: " << data.recordCount() << "\n";
    return 0;
}