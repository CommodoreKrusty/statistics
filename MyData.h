//Generally templates don't have a source file so everything is in the headder. Might change that in the future.
#include <stdarg.h>
#include <list>
#include <array>
#include <iostream>
#include <string>
#include <cctype>
#include <cassert>

//using namespace std;
 using std::string,
	std::size_t,
	std::array,
	std::list,
	std::initializer_list,
	std::copy,
	std::to_string,
	std::cerr;

/* template <typename t, size_t s> class myArray: public array<t, s>{
    public:
    myArray(){};
    myArray(initializer_list<t> l){

        assert(l.size() == s);
        copy(l.begin(), l.end(), this->begin());
    }

    void setColumnToSort(int col){
        this->col_to_sort = col;
    };

    bool operator<(myArray<t, s> a){
        if (this->at(col_to_sort) < a[col_to_sort]){
            return true;
        }
        return false;
    };

    int col_to_sort;
};


template <class T, int size> class newList : public list<T>{
	public:
		void mysort(int col);
//		array<newList<T, size>, 2> splitList(int col);
};


template <class T, int size> void newList<T, size>::mysort(int col){
	typename newList<T, size>::iterator i = this->begin();
	int c = 0;
	while (i != this->end()){
		i->setColumnToSort(col);
		++i;
	}

	this->sort();
} 
*/
template <typename t, size_t s> class myArray: public array<t, s>{
	public:
	myArray(){};
	myArray(initializer_list<t> l){

		assert(l.size() == s);
		copy(l.begin(), l.end(), this->begin());
	};

	void setColumnToSort(int col){
		this->col_to_sort = col;
	};

	bool operator<(myArray a){
		if (this->at(col_to_sort) < a[col_to_sort]){
			return true;
		}
		return false;
	};

	int col_to_sort;
};


template <typename t> class newList : public list<t>{
	public:
		void mysort(int col){
			typename list<t>::iterator i = this->begin();
			int c = 0;
			while (i != this->end()){
				i->setColumnToSort(col);
				++i;
			}
		
			this->sort();
		};
//		array<newList<T, size>, 2> splitList(int col);
};

#define LOCAL
#ifdef LOCAL

template <typename t, size_t s> class MyData{
//	private:
//		int sort_col;
	public:

		
		/*template <class T, int size> array<newList<T, size>, 2> newList<T, size>::splitList(int col){
			array<newList<T, size>, 2> a;
			newList<T, size> l1, l2;
			return a;
		}
		*/
		
/* 		template <t, size_t> void newList<t, size_t>::mysort(int col){
			typename newList<t, size>::iterator i = this->begin();
			int c = 0;
			while (i != this->end()){
				i->setColumnToSort(col);
				++i;
			}
		
			this->sort();
		}
 */	
		size_t record_size = s;//number of columns
		typedef myArray<t, s> newArray;
		typedef newList<newArray> dataList;
//		typedef array<t,s> newArray;
//		typedef list<newArray> dataList;
		dataList data;
//		static bool compare_val(myArray first, myArray second);


	public:

		typedef typename dataList::iterator iterator;
//		typedef datalist::const_iterator const_iterator;
		iterator begin(){return data.begin();}
		iterator end(){return data.end();}
		newArray front(){return data.front();}
		newArray back(){return data.back();}

		MyData(){};
		~MyData(){
			while(!data.empty()){
			//	don't know if I need to delete members of the array here. I don't think so.
			//	myarray a = data.front();
				data.pop_front();
			}
				
		};
//		void addRecord(int c, ...);// 
		void addRecord(newArray a){
			data.push_back(a);
		};
		void print(){
			string msg;
			iterator i = begin();
			int c = 0;
			while (i != end()){
				myArray a = *i;
				for(typename array<t, s>::iterator b = a.begin(); b != a.end(); ++b){
					msg += to_string(*b);
					msg += "\t";
				}
				msg += '\n';
				++c;
				++i;
			}
			cerr << msg;
		
		};

		int recordCount(){
			return data.size();
		};

		void sort(int col){
			try{
				if (col >= 0 and col < s){
					data.mysort(col);
				}
			}catch(int error){
				cerr << "out of bounds" << "\n";//exceptionHandler(error);
			}		
		};
		

//		friend ostream& operator<<(ostream &s, MyData n)<T>;


};

/* template <class T, int size> MyData<T, size>::MyData(){
}
 */

/* template <class T, int size> MyData<T, size>::~MyData()
{
	while(!data.empty()){
//		don't know if I need to delete members of the array here. I don't think so.
//		myarray a = data.front();
		data.pop_front();
	}
}
 */
/* template <class T, int size> void MyData<T, size>::addRecord(newArray a){
	data.push_back(a);
}
 */

//variadic function
/*template <typename T, int size> void MyData<T, size>::addRecord(int c, ...){
	myArray a;
	int result = 0;
	va_list args;
	va_start(args, c);
	for (int i = 0; i < record_size; ++i) {
		result = va_arg(args, T);
		a[i] = result;
	}
	va_end(args);
	addRecord(a);
}
*/
/* template <class T, int size> int MyData<T, size>::recordCount(){
	return data.size();
}
 */

/* template <class T, int size> void MyData<T, size>::print(){
//	cout << "size: " << to_string(recordCount()) << endl;
	string msg;
	iterator i = begin();
	int c = 0;
	while (i != end()){
		myArray a = *i;
		for(typename array<T, size>::iterator b = a.begin(); b != a.end(); ++b){
			msg += to_string(*b);
			msg += "\t";
		}
		msg += '\n';
		++c;
		++i;
	}
	cerr << msg;
} */

/*template <typename T, int size> bool MyData<T, size>::compare_val(myArray first, myArray second){
	return false;
}*/

/* template <class T, int size> void MyData<T, size>::sort(int col){
//	sort_col = col;
	try{
		if (col >= 0 and col < size){
			data.mysort(col);
		}
	}catch(int error){
		cerr << "out of bounds" << "\n";//exceptionHandler(error);
	}

} */

/*template <typename T> ostream& operator<<(ostream &s, MyData n){
	string msg;
	MyData::iterator it;
	for (it = begin(); it != end(); ++it){
		msg += to_string(*it);
		msg += "*";
	}

	msg.pop_back();
	return s << msg;

}*/

#endif
