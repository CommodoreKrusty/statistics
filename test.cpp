#include <iostream>
#include <random>
#include <time.h>
#include <list>
#include <algorithm>
#include <cstdlib>
#include "/home/paul/Documents/c++/statistics/MyStatistics.h"

using std::mt19937, std::div, std::div_t, std::list, std::for_each, std::advance, std::cerr;

template < class Key, int T> class SuperMyStatistics : public MyStatistics< Key, T>{
    public:
    SuperMyStatistics(){}
    ~SuperMyStatistics(){
//        this->data.empty();
    }

};

typedef SuperMyStatistics<int, 2> dataset;
typedef list<dataset> datasets;

void FillDataSet(mt19937 rnd, dataset* s);

int main(){
    mt19937 rnd(time(NULL));
    datasets mydata;
    dataset s;
    FillDataSet(rnd, &s);
    mydata.push_back(s);

    auto print = [](dataset v){
        v.print();
    };
    s.multisort({0, 1});
    print(s);
//    for_each(mydata.begin(), mydata.end(), print);

/*    datasets::iterator i = mydata.begin();
    while(i != mydata.end()){
        i->sort();
        advance(i, 1);
    }
*/
//    for_each(mydata.begin(), mydata.end(), print);

    return 0;
}

void FillDataSet(mt19937 rnd, dataset* s){
    for(int i = 0; i < 50; i++){
        dataset::MyArray j;

//        div_t a = div(rnd() % 100 + 1, 10);
//        div_t b = div(rnd() % 100 + 1, 10);
        unsigned long a = rnd() % 10 + 1;
        unsigned long b = rnd() % 10 + 1;

//        j[0] = (double) a.quot + (double)a.rem / 10;
//        j[1] = (double) b.quot + (double)b.rem / 10;
        j[0] = a;
        j[1] = b;

        s->addRecord(j);
    }

}
