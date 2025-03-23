#include <iostream>
#include <list>
#include <vector>
#include <algorithm>

using std::list, std::vector, std::for_each, std::advance, std::cerr;

template < class Key> class SuperList{
public:
    template <class S> class sortvector : public vector <S>{
    };

    typedef sortvector<int> myvector;

    myvector data;

    typename myvector::iterator begin(){return data.begin();}
    typename myvector::iterator end(){return data.end();}

    SuperList(){}
    ~SuperList(){

    };

    void add(Key n){
        data.push_back(n);
    };

};

template < class Key> class SuperDuperList: public SuperList<Key>{
public:
    SuperDuperList(){}
    void print(){
/*        typename SuperDuperList::myvector::iterator i = this->begin();
        while(i != this->end()){
            cerr << *i << " ";
            advance(i, 1);
        }
*/
        for_each(this->begin(), this->end(), [](int i){ cerr << i << " ";});
        cerr << "\n";
    };
};

typedef SuperDuperList<int> mylist;

int main(){
    mylist l;
    for(int x = 0; x < 50; x++){
        l.add(x);
    }
    l.print();
    return 0;
}
