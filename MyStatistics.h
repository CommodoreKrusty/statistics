#include <stdarg.h>
#include <map>
#include <cmath>
#include <list>
#include <iostream>
#include <array>
#include <string>
#include "MyData.h"

//using namespace std;
using std::string,
	std::array,
	std::map,
	std::list,
	std::to_string,
	std::cerr;

//don't screw with these
#define ERROR_OUT_OF_RANGE 1

//Default values. If you want to change any of these values do it here and not down in the code.
#define DEGREES_OF_FREEDOM 1
#define ALPHA_LEVEL 0.05
#define ONE_TAILED 0
#define TWO_TAILED 1

template < class Key, class T> class newMap : public map< Key, T>{
};

class RegressionLineEquation{
	double mx, bx;
	public:
		RegressionLineEquation(double m, double b);
		string getEquation();
		double getYIntercept(){return bx;}
		double getSlope(){return mx;}
};

RegressionLineEquation::RegressionLineEquation(double m, double b){
	mx = m;
	bx = b;
}

string RegressionLineEquation::getEquation(){
	string msg;
	if (bx < 0){
		msg = "y = " + to_string(mx) + "x - " + to_string(bx*(-1));
	}else if(bx > 0){
		msg = "y = " + to_string(mx) + "x + " + to_string(bx);
	}else{
		msg = "y = "+ to_string(mx) + "x"; 
	}
	return msg;
}

/*
* 'size' is the number of fields in a record.
*/
template <typename T, int size> class MyStatistics: public MyData<T, size>{

	public:
		typedef array<T, size> myArray;
		typedef newList<myArray> dataList;
		typedef typename dataList::iterator iterator;
		iterator begin(){return this->begin();}
		iterator end(){return this->end();}

		MyStatistics();
		~MyStatistics();
		int exceptionHandler(int error);
		myArray totalColumns(); 
		T totalColumn(int col);
		double colMean(int col);
		T colMedian(int col);
		T colRange(int col);
		T colMidrange(int col);
		T IQR(int col);
		double variance(int col, int sample = 1);//check this
		double sandarddev(int col, int sample = 1);
		double percentile(int col, T num);
		double MAD(int col);
		map<T, int> colMode(int col);
		double meanofSquared(int col);
		double CorrelationCoefficient(int col1, int col2);
		RegressionLineEquation LineOfRegression(int col1, int col2);
		list<array<double,2>> residuals(int col1, int col2);
		RegressionLineEquation squaredErrorofLineofRegression(int colX, int colY);
		void printfrequencyTable(int col);
		newMap<T, int> frequencyTable(int col, bool group = false);
		double RSquared(int colX, int colY);

		bool zScore(T num, int col = 0, int sample = 1, double alpha = ALPHA_LEVEL);
//		bool zScoreTest(T num, double colMean, double stdev, int sample_size, double alpha = ALPHA_LEVEL);
		bool zScoreTest(T Ha, T Ho, int sample_size, double alpha = ALPHA_LEVEL);//, double stdev
		bool zScoreTest(T Ha, T Ho, double stdev, double alpha = ALPHA_LEVEL);//, double stdev
		double zProbabilityDensity(double x);
		double ztrapezoid(double x, double s);
		double zScorePValue(double z);

		double zStatistic(T num, int col = 0, int sample = 1);
		double zTotScore(double z);
//		double pValue(int zStat = 0);

		bool tTest(int col = 0, double Ha = 0, int tailed = ONE_TAILED, double alpha = ALPHA_LEVEL);
		double tWhereToStartCounting(int df);
		double ProbabilityDistributionFormula(double t, int df = DEGREES_OF_FREEDOM);


//		void tValueTest(double t, double alpha = ALPHA_LEVEL);
		bool tHypothesisTest(int col = 0, double Ha = 0, int tailed = ONE_TAILED, double alpha = ALPHA_LEVEL);
		double tStatistic(int col, double Ha = 0);
		bool tTable(double t, int tailed = ONE_TAILED, double alpha = ALPHA_LEVEL);


		bool chiSquaredTest(int colExpeced = 0, int colActual = 1, double alpha = ALPHA_LEVEL);
		

		protected:
		double meanXY(int colX, int colY);
		double meanOfColumnSquared(int col);

		double TotalSEwithLine(int colX, int colY);
//		double SEwithLine();
		double TotalSEfromMean(int col);
		double trapezoid(double x = 0.0, double s = 1.0);


//		double ProbabilityDistributionFormula(double t, int df = DEGREES_OF_FREEDOM);
		double tValue(double tStat);
		double Ttrapezoid(double t, double dx, int df = DEGREES_OF_FREEDOM);
//		double tStatistic(int col, double populationMean = 0);
		double tStatistic(double Ho, double Ha, double sampleStandardDeviation, int sampleSize);//work on this stuff
		double Ctrapezoid(double t, double dx, int df = DEGREES_OF_FREEDOM);

//chiSquared stuff
		double chiSquared(int colExpeced = 0, int colActual = 1);
		double chiSquaredIntegral(double cX = 1, double df = 0, double x = 0);
		double chiSquaredProbabilityDensity(double df, double x = 0);


	private:
		void areaForAllCurves();
		double halfTotalAreaUnderCurve(int df = DEGREES_OF_FREEDOM);
		double halfAreaUnderTDistribution(int df = DEGREES_OF_FREEDOM);

};
//notice it's .tpp and not .cpp just to remind me it's template stuff.
//#include "MyStatistics.tpp"//this is a strange place to include the source file but it doesn't work if it's included at the top

template <class T, int size> MyStatistics<T, size>::MyStatistics(){
//	this->addRecord({97,36,43,22,76});//I don't know why 'this->' is necessary. It wouldn't be if it wasn't a template.
}

template <class T, int size> MyStatistics<T, size>::~MyStatistics()
{
}

template <class T, int size> array<T, size> MyStatistics<T, size>::totalColumns(){
	array<T, size> a;

	for(int n = 0; n < size; ++n){
		a[n] = (T)totalColumn(n);
	}

	return a;
}

template <class T, int size> T MyStatistics<T, size>::totalColumn(int col){
	T t = (T)0;
	try{
		if (col >= 0 and col < size){
			typename MyData<T, size>::iterator i = MyData<T, size>::begin();
			while (i != MyData<T, size>::end()){
				myArray m = *i;
				t += m.at(col);
				++i;
			}
		}
		else{
			throw ERROR_OUT_OF_RANGE;//return -1;//error. Maybe throw exception?
		}
	}
	catch(int error){
		exceptionHandler(error);
	}
	return t;
}

template <class T, int size> double MyStatistics<T, size>::colMean(int col){
	double a = 0.0;
	try{
		if (col >= 0 and col < size){
			double total = totalColumn(col);
			a = total / (double)this->recordCount();
		}else{
			throw ERROR_OUT_OF_RANGE;//return -1;//error. Maybe throw exception?
		}
	}catch(int error){
		exceptionHandler(error);
	}
	return a;
}

template <class T, int size> double MyStatistics<T, size>::meanofSquared(int col){
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	double t = 0.0;
	while (i != MyData<T, size>::end()){
		myArray m = *i;
		double n = (double)m.at(col);
		double p = pow(n, 2);
		t += p;
		++i;
	}

	return t / (double)this->recordCount();
}
/*
template <class T, int size> double MyStatistics<T, size>::meanofSquared(int col){
	typename list<myArray>::iterator i = this->begin();
	double t = 0.0;
	while (i != this->end()){
		myArray m = *i;
		double n = (double)m.at(col);
		double p = pow(n, 2);
		t += p;
		++i;
	}

	return t / (double)this->recordCount();
}*/


template <class T, int size> map<T, int> MyStatistics<T, size>::colMode(int col){
	map<T, int> m;
	int max_val = 0;
	map<T, int> section = frequencyTable(col);
	typename map<T, int>::iterator i = section.begin();
//	max_val = i->second;
	while (i != section.end()){
//		cerr << i->first << " => " << i->second << '\n';
		if (i->second > max_val){
			max_val = i->second;
		}
		++i;
	}
	i = section.begin();
	while (i != section.end()){
		if (i->second == max_val){
			m[i->first] = i->second;
		}
		++i;
	}
	return m;
}

template <class T, int size> void MyStatistics<T, size>::printfrequencyTable(int col){
	cerr << "col: " << col << '\n';
	map<T, int> section = frequencyTable(col);
	typename map<T, int>::iterator i = section.begin();
	while (i != section.end()){
		cerr << i->first << " => " << i->second << '\n';
		++i;
	}
}

template <class T, int size> newMap<T, int> MyStatistics<T, size>::frequencyTable(int col, bool group){
	newMap<T, int> section;
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		myArray m = *i;
		T t = m[col];
		typename map<T, int>::iterator si = section.find(t);
		if (si != section.end()){
			section[t] += 1;
			
		}else{
			section[t] = 1;
		}
		++i;
	}
	return section;
}

//not entirely certain how to calculate median.
template <class T, int size> T MyStatistics<T, size>::colMedian(int col){
	T median = 0;
	try{
		if (col >= 0 and col < size){
			MyData<T, size>::sort(col);
			int l = MyData<T, size>::recordCount() - 1;
			int p = l%2;
			cerr << "p: " << to_string(p) << '\n';
			int n = 0;
			if (p == 0){
				n = (int)l/2;
			}else{
				n = (l+1)/2;
			}
			typename MyData<T, size>::iterator i = MyData<T, size>::begin();
			int c = 0;
			while (c < n && i != MyData<T, size>::end()){
				++c;
				++i;
			}
			myArray t = *i;
		 	median = t[col];
		}else{
			throw ERROR_OUT_OF_RANGE;//return -1;//error. Maybe throw exception?
		}
	}catch(int error){
		exceptionHandler(error);
	}


	return median;
}

//The highest value in a column minus the lowest value
template <class T, int size> T MyStatistics<T, size>::colRange(int col){
	T mode = (T)0;
	MyData<T, size>::sort(col);
	myArray t = MyData<T, size>::front();
	T v = t[col];
	t = MyData<T, size>::back();
	T w = t[col];
	return (T)(w - v);
}


//The average of the lowest and highest elements in a column
template <class T, int size> T MyStatistics<T, size>::colMidrange(int col){
	T r = colRange(col);
	return r/(T)2;
}

//Interquartile Range
template <class T, int size> T MyStatistics<T, size>::IQR(int col){
	MyData<T, size>::sort(col);
	T Q1, Q2;
	int pos1, pos2;

	int count = MyData<T, size>::recordCount();
	pos1 = (count/2)/2;
	pos2 = (count/2)+(count/4);
	
//	cerr << "count: " << to_string(count)<< " pos1: " << to_string(pos1) << " pos2: " << to_string(pos2) << '\n';

	int c = 0;
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		if(c == pos1){
			myArray t = *i;
		 	Q1 = t[col];
		}else if (c == pos2){
			myArray t = *i;
		 	Q2 = t[col];
		}
		++i;
		++c;
	}

//	cerr << "Q2: " << to_string(Q2) << " Q1: " << to_string(Q1) << '\n';

	return (T)(Q2-Q1);
}

//standard deviation
template <class T, int size> double MyStatistics<T, size>::sandarddev(int col, int sample){
	return sqrt(variance(col, sample));
}

//'sample' defaults to 1 for sample.
template <class T, int size> double MyStatistics<T, size>::variance(int col, int sample){
	if (sample != 1){
		sample = 0;//Anything other than 1 will become 0 for population
	}
	double a = 0;
	double m = colMean(col);
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		myArray l = *i;
		a += pow((double)l[col] - m, 2);
		++i;
	}
	
	int count = MyData<T, size>::recordCount() - sample;
	a /= count;
	return a;
}

template <class T, int size> double MyStatistics<T, size>::percentile(int col, T num){
	MyData<T, size>::sort(col);
	int count = 0;
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		myArray l = *i;
//		T t = l[col];
		if (l[col] < num){
			count += 1;
		}
		++i;
	}
	cerr << "count: " << to_string(count) << '\n';
	int rc = MyData<T, size>::recordCount();
	cerr << "recordCount: " << to_string(rc) << '\n';
	double p = (double)count / (double)rc;
	return p;
}

template <class T, int size> double MyStatistics<T, size>::MAD(int col){
	double total = 0;
	double m = colMean(col);
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		myArray l = *i;
		T t = l[col] - m;
		if (t < 0){
			t *= -1;
		}
		total += (double) t;
		++i;
	}
	double t = (double)total/(double)MyData<T, size>::recordCount();
	return t;
}

template <class T, int size> double MyStatistics<T, size>::CorrelationCoefficient(int col1, int col2){
	double mean1 = colMean(col1);
	double stdev1 = sandarddev(col1);
	double mean2 = colMean(col2);
	double stdev2 = sandarddev(col2);

	double coef = (double)1/((double)MyData<T, size>::recordCount() - 1.0);
	int temp = 0;
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		myArray l = *i;
		T a = l[col1] - mean1;
		T b = l[col2] - mean2;
		double t = (double)(a * b);
		if (t != 0){
			temp += 1;
		}
		++i;
	}

	return coef * (double)temp / (stdev1 * stdev2);
}

template <class T, int size> RegressionLineEquation MyStatistics<T, size>::LineOfRegression(int col1, int col2){
	double mean[2];
	mean[0] = colMean(col1);
	mean[1] = colMean(col2);

	double stdev[2];
	stdev[0] = sandarddev(col1);
	stdev[1] = sandarddev(col2);

	double r = CorrelationCoefficient(col1, col2);

	double m = r*stdev[1]/stdev[0];
	double b = mean[1] - r*stdev[1]/stdev[0]*mean[0];

	return RegressionLineEquation(m, b);
}

template <class T, int size> list<array<double,2>> MyStatistics<T, size>::residuals(int col1, int col2){
	list<array<double,2>> l;
	RegressionLineEquation r = LineOfRegression();
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		myArray temp = *i;
		array<double, 2> a = {(double)temp[col1], (double)(temp[col2] - (r.getSlope()*temp[col2]) + r.getYIntercept())};
		l.push_back(a);
		++i;
	}

	return l;
}
//Squared Error of the Line of Regression
template <class T, int size> RegressionLineEquation MyStatistics<T, size>::squaredErrorofLineofRegression(int colX, int colY){
	double meanX = colMean(colX);
	double meanY = colMean(colY);
	double mXY = meanXY(colX, colY);
	double mXsq = meanOfColumnSquared(colX);
	double mXmY = meanX * meanY;
	double numerator = mXmY - mXY;
	double denominator = pow(meanX, 2) - mXsq;
	double m = numerator/denominator;
	double t = m * meanX;
	double b = meanY - t;

	return RegressionLineEquation(m,b);
}

template <class T, int size> double MyStatistics<T, size>::meanXY(int colX, int colY){
	double total = 0.0;
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		myArray temp = *i;
		total += temp[colX] * temp[colY];
		++i;
	}
	return total/(double)MyData<T, size>::recordCount();
}

template <class T, int size> double MyStatistics<T, size>::meanOfColumnSquared(int col){
	double total = 0.0;
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		myArray temp = *i;
		total += pow(temp[col],2);
		++i;
	}
	return total / (double)MyData<T, size>::recordCount();
}
//END Squared Error

//R-Squared stuff
template <class T, int size> double MyStatistics<T, size>::RSquared(int colX, int colY){
	double tSE = TotalSEwithLine(colX, colY);
	double tSEm = TotalSEfromMean(colY);
//	Percent Of Total Variation Not Explained by Variation In X
	double pTVNot = tSE/tSEm;
	return 1.0 - pTVNot;
}

template <class T, int size> double MyStatistics<T, size>::TotalSEwithLine(int colX, int colY){
	double total = 0.0;
	RegressionLineEquation r = squaredErrorofLineofRegression(colX, colY);
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		myArray temp = *i;
		//actual Y
		T y = temp[colY];
		//expected Y
		double eY = r.getSlope() * temp[colX] + r.getYIntercept();
		double differenceY = (double)y - eY;
		differenceY = pow(differenceY, 2.0);
		total += differenceY;
		++i;
	}

	return total;
}

template <class T, int size> double MyStatistics<T, size>::TotalSEfromMean(int col){
	double total = 0.0;
	double meanCol = colMean(col);
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	while (i != MyData<T, size>::end()){
		myArray temp = *i;
		T y = temp[col];
		double t = (double)y - meanCol;
		total += pow(t, 2.0);
		++i;
	}
	return total;
}
//End of R-Squared stuff

//Z score stuff

template <class T, int size> bool MyStatistics<T, size>::zScoreTest(T Ha, T Ho, double stdev, double alpha){//, double stdev
	cerr << "here Ha:" << to_string(Ha) << " Ho: " << to_string(Ho) << '\n';
	double z = (Ha - Ho)/stdev;//(stdev/sqrt(sample_size));
	if (z < 0){
		z *= -1;
	}

	cerr << "z: " << to_string(z) << '\n';
	
	double area = zScorePValue(z);

	if (area < alpha){
		return false;//reject null hypothesis
	}else{
		return true;//fail to reject null hypothesis
	}
}


template <class T, int size> bool MyStatistics<T, size>::zScoreTest(T Ha, T Ho, int sample_size, double alpha){//, double stdev
	double z = ((double)Ha/sample_size - (double)Ho)/sqrt((Ho*(1-Ho))/sample_size);//(stdev/sqrt(sample_size));
	if (z < 0){
		z *= -1;
	}

	cerr << "z: " << to_string(z) << '\n';

	double area = zScorePValue(z);

	if (area < alpha){
		return false;//reject null hypothesis
	}else{
		return true;//fail to reject null hypothesis
	}
}

template <class T, int size> bool MyStatistics<T, size>::zScore(T num, int col, int sample, double alpha){
	if (sample != 1){
		sample = 0;//Anything other than 1 will become 0 for population
	}
	double meanCol = colMean(col);
	double stdev = sandarddev(col, sample);
	return zScoreTest(num, meanCol, stdev, alpha);
/*	if (sample != 1){
		sample = 0;//Anything other than 1 will become 0 for population
	}

	double meanCol = colMean(col);
	double stdev = sandarddev(col, sample);

	double z = ((double)num - (double)meanCol)/stdev;
	if (z < 0){
		z *= -1;
	}
	cerr << "z: " << to_string(z) << '\n';
	
	double area = 0;
	double dx = (double)1/pow(10, 4);
	double x = 38.575520;//0;//

	do{
		area += ztrapezoid(x, dx);
		x -= dx;
	}while(x >= z);

//	cerr << "area: " << to_string(area) << '\n';
	if (area < alpha){
		return false;//reject null hypothesis
	}else{
		return true;//fail to reject null hypothesis
	}
*/
}

template <class T, int size> double MyStatistics<T, size>::zScorePValue(double z){
	double area = 0;
	double dx = (double)1/pow(10, 4);
	double x = 38.575520;

	do{
		area += ztrapezoid(x, dx);
		x -= dx;
	}while(x >= z);
	return area;
}


template <class T, int size> double MyStatistics<T, size>::zProbabilityDensity(double x){
	const double pi = 3.14159265358979323846;
	const double E = 2.71828182845904523536;

	double a = pow(E,(-1.0)*(pow(x,2.0))/2.0)/pow(2.0 * pi, 0.5);

	return a;
}

template <class T, int size> double MyStatistics<T, size>::ztrapezoid(double x, double s){
	double a = zProbabilityDensity(x);
	double b = zProbabilityDensity(x - s);
	return (a + b)/2 * s;
}

#ifdef DONTSKIP
template <class T, int size> double MyStatistics<T, size>::zStatistic(T num, int col, int sample){
	if (sample != 1){
		sample = 0;//Anything other than 1 will become 0 for population
	}

	double meanCol = colMean(col);
	double stdev = sandarddev(col, sample);

	return ((double)num - (double)meanCol)/stdev;
}

/*template <class T, int size> double MyStatistics<T, size>::pValue(int zStat){
	double area = 0;
	double dx = (double)1/20000;
	double n = 0.0;
	while((n + dx)<= fabs(zStat)){
		area += trapezoid( n, dx);
		n += dx;
	}
	if (zStat >= 0){
		return (area + 0.5);
	}

	return (0.5 - area);
}*/

template <class T, int size> double MyStatistics<T, size>::trapezoid(double x, double s){
	const double pi = 3.14159265358979323846;
	const double E = 2.71828182845904523536;

	double a = pow(E,(-1.0)*(pow((double)x,2.0))/2.0)/pow(2.0 * pi, 0.5);
	double b = pow(E, (-1.0)*(pow((double)x + (double)s, 2.0)/2.0))/pow(2.0 * pi, 0.5);
	double area = ((a + b)/2)*s;
 	return area;
}

#endif

//Found this on the web. does this even work?
template <class T, int size> double MyStatistics<T, size>::zTotScore(double z){
	return (z * 10.0) + 50.0;
}

//End Z score stuff

//t statistic stuff for sample size < 30

template <class T, int size> bool MyStatistics<T, size>::tTest(int col, double Ha, int tailed, double alpha){
	double df = MyData<T, size>::recordCount() - 1;
	double t = tStatistic(col, Ha);
	double tVal = tValue(t);
//	cerr << "t " << to_string(t) << '\n';
//	cerr << "tVal: " << to_string(tVal) << '\n';
	if (tVal < alpha){
		cerr << "reject Ho" << '\n';
		return false;//reject null hypothesis.
	}else{
		cerr << "fail to reject Ho" << '\n';
		return true;//accept null hypothesis.
	}

/*	switch(tailed)
	{
		case ONE_TAILED:
		{
			if (tVal < alpha){
				cerr << "reject Ho" << '\n';
				return false;//reject null hypothesis.
			}else{
				cerr << "fail to reject Ho" << '\n';
				return true;//accept null hypothesis.
			}
		}
		case TWO_TAILED:
		{
		}
	}
*/
//	cerr << "area: " << to_string(area) << '\n';

	return false;
/*
	double df = MyData<T, size>::recordCount() - 1;
	double t = tStatistic(col, populationMean);
	double tVal = tValue(t);
	double half_area = 0.5;//halfAreaUnderTDistribution(df);
	switch(tailed)
	{
		case ONE_TAILED:
		{
			double j = tVal + half_area;
			double b = 2.0 * half_area;
			double z = (b - j)/b;

			if (z < alpha){
				return false;//reject null hypothesis.
			}else{
				return true;//accept null hypothesis.
			}

			return true;
		}
		case TWO_TAILED:
		{
//		double t = tStatistic(col, populationMean);
//		double tVal = tValue(t);
//		double df = MyData<T, size>::recordCount() - 1;
			double bound = half_area * (1 - (alpha/2));

			if (tVal > bound){//(upper > (1.0 - alpha)){
				cerr << " reject the NULL hypothesis" << '\n';
				return false;
			}
			cerr << " accept" << '\n';
			return true;
			break;
		}
	}
	return true;
*/
}

template <class T, int size> bool MyStatistics<T, size>::tHypothesisTest(int col, double Ha, int tailed, double alpha){
	double t = tStatistic(col, Ha);
	bool r = tTable(t, tailed, alpha);
	if(r){
		cerr << "fail to reject Ho" << '\n';
	}else{
		cerr << "reject Ho" << '\n';
	}
	return r;
}

template <class T, int size> double MyStatistics<T, size>::tStatistic(int col, double Ha){
	double meanSample = colMean(col);
//	cerr << "meanSample: " << to_string(meanSample) << '\n';

//	double SE = meanSample/sqrt(count);
//	double T = 0;
//	double meanPopulation = meanSample + T*SE; 
	double stdevSample = sandarddev(col);
//	cerr << "stdevSample: " << to_string(stdevSample) << '\n';
	int count = MyData<T, size>::recordCount();
//	double t = (double)(meanSample - populationMean)/(stdevSample/sqrt(count));
	
	double t = tStatistic(meanSample, Ha, stdevSample, count);
	cerr << "t: " << to_string(t) << '\n';
	return t;
}

template <class T, int size> double MyStatistics<T, size>::tStatistic(double Ho, double Ha, double sampleStandardDeviation, int sampleSize){
	return (Ho - Ha)/(sampleStandardDeviation/sqrt(sampleSize));
}
/*
template <class T, int size> double MyStatistics<T, size>::tStatistic(double sampleMean, double populationMean, double sampleStandardDeviation, int sampleSize){
	return (sampleMean - populationMean)/(sampleStandardDeviation/sqrt(sampleSize));
}*/

/*template <class T, int size> void MyStatistics<T, size>::tValueTest(double t, double alpha){
	double df = 1;//MyData<T, size>::recordCount() - 1;
	double area = 0.0;
	double dx = (double)1/pow(10, 5);
	double x = tWhereToStartCounting(df);// + 120.0
	while(area <= 0.05){
		area += Ttrapezoid( x, dx, df);
//		cerr << "area: " << to_string(area) << '\n';
/*		if(area == 0.05){
			cerr << "t score: " << to_string(x) << '\n';
		}*/
/*		x -= dx;
	}
	cerr << "t score: " << to_string(x) << '\n';
}*/

/*
* tStat is output from tStatistic();
*/
template <class T, int size> double MyStatistics<T, size>::tValue(double tStat){
	double df = MyData<T, size>::recordCount() - 1;
	double area = 0.0;
	double dx = (double)1/pow(10, 4);//if you want dx to be smaller make 4 larger, but it will be slower.
	double x = tWhereToStartCounting(df);
	while(x > tStat){
		area += Ttrapezoid( x, dx, df);
		x -= dx;
	}
	return area;
}

/*
template <class T, int size> double MyStatistics<T, size>::tValue(double tStat){
	double df = MyData<T, size>::recordCount() - 1;
	double area = 0;
	double dx = (double)1/pow(10, 4);
	double n = 0.0;
	while( (n + dx) <= fabs(tStat)){//(n + dx)
		area += Ttrapezoid( n, dx, df);
		n += dx;
	}
	return area;
}
*/
template <class T, int size> double MyStatistics<T, size>::Ttrapezoid(double t, double dx, int df){
	double a = ProbabilityDistributionFormula(t, df);
	double b = ProbabilityDistributionFormula(t - dx, df);
	return ((a + b)/2 * dx);
}

//https://www.thoughtco.com/students-t-distribution-formula-3126276
//Student's t Distribution Formula:
//v - degrees of freedom
//gamma((v + 1)/2)/(sqrt(v*pi)*gamma(v/2))*(1+t^2/v)^(-(v+1)/2)
template <class T, int size> double MyStatistics<T, size>::ProbabilityDistributionFormula(double t, int df){
	const double pi = 3.14159265358979323846;
	//v - degrees of freedom
	double v = (double)df;//MyData<T, size>::recordCount() - 1;//
	double g = tgamma((v + 1.0)/2.0);
	double d = sqrt(v * pi)*tgamma(v/2.0);
	double c = 1.0 + pow(t,2)/v;
	double p = -((v+1.0)/2.0);
	return (g/d)*pow(c,p);
}


//I think I'm overthinking this.
template <class T, int size> bool MyStatistics<T, size>::tTable(double t, int tailed, double alpha){
//#define ONE_TAILED 0
//#define TWO_TAILED 1
	int df = MyData<T, size>::recordCount() - 1;

	if (alpha == 0.05){
		switch(df){
		case 1:
			if(tailed == ONE_TAILED){
				if (t > 6.313751){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -12.7062 || t > 12.7062){
					return false;
				}
			}else{
			}
		case 2:
			if(tailed == ONE_TAILED){
				if (t > 2.919986){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -4.302653 || t > 4.302653){
					return false;
				}
			}else{
			}
		case 3:
			if(tailed == ONE_TAILED){
				if (t > 2.353363){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -3.182446 || t > 3.182446){
					return false;
				}
			}else{
			}
		case 4:
			if(tailed == ONE_TAILED){
				if (t > 2.131802){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.776298 || t > 2.776298){
					return false;
				}
			}else{
			}
		case 5:
			if(tailed == ONE_TAILED){
				if (t > 2.015036){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.570543 || t > 2.570543){
					return false;
				}
			}else{
			}
		case 6:
			if(tailed == ONE_TAILED){
				if (t > 1.943176){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.446899 || t > 2.446899){
					return false;
				}
			}else{
			}
		case 7:
			if(tailed == ONE_TAILED){
				if (t > 1.894577){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.364619 || t > 2.364619){
					return false;
				}
			}else{
			}
		case 8:
			if(tailed == ONE_TAILED){
				if (t > 1.859547){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.306002 || t > 2.306002){
					return false;
				}
			}else{
			}
		case 9:
			if(tailed == ONE_TAILED){
				if (t > 1.833113){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.262156 || t > 2.262156){
					return false;
				}
			}else{
			}
		case 10:
			if(tailed == ONE_TAILED){
				if (t > 1.812461){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.228138 || t > 2.228138){
					return false;
				}
			}else{
			}
		case 11:
			if(tailed == ONE_TAILED){
				if (t > 1.795885){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.200985 || t > 2.200985){
					return false;
				}
			}else{
			}
		case 12:
			if(tailed == ONE_TAILED){
				if (t > 1.782288){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.178813 || t > 2.178813){
					return false;
				}
			}else{
			}
		case 13:
			if(tailed == ONE_TAILED){
				if (t > 1.770933){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.160369 || t > 2.160369){
					return false;
				}
			}else{
			}
		case 14:
			if(tailed == ONE_TAILED){
				if (t > 1.76131){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.144787 || t > 2.144787){
					return false;
				}
			}else{
			}
		case 15:
			if(tailed == ONE_TAILED){
				if (t > 1.75305){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.131449 || t > 2.131449){
					return false;
				}
			}else{
			}
		case 16:
			if(tailed == ONE_TAILED){
				if (t > 1.745884){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.119905 || t > 2.119905){
					return false;
				}
			}else{
			}
		case 17:
			if(tailed == ONE_TAILED){
				if (t > 1.739607){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.109816 || t > 2.109816){
					return false;
				}
			}else{
			}
		case 18:
			if(tailed == ONE_TAILED){
				if (t > 1.734064){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.100922 || t > 2.100922){
					return false;
				}
			}else{
			}
		case 19:
			if(tailed == ONE_TAILED){
				if (t > 1.729133){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.093024 || t > 2.093024){
					return false;
				}
			}else{
			}
		case 20:
			if(tailed == ONE_TAILED){
				if (t > 1.724718){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.085963 || t > 2.085963){
					return false;
				}
			}else{
			}
		case 21:
			if(tailed == ONE_TAILED){
				if (t > 1.720743){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.079614 || t > 2.079614){
					return false;
				}
			}else{
			}
		case 22:
			if(tailed == ONE_TAILED){
				if (t > 1.717144){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.073873 || t > 2.073873){
					return false;
				}
			}else{
			}
		case 23:
			if(tailed == ONE_TAILED){
				if (t > 1.713872){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.068658 || t > 2.068658){
					return false;
				}
			}else{
			}
		case 24:
			if(tailed == ONE_TAILED){
				if (t > 1.710882){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.063899 || t > 2.063899){
					return false;
				}
			}else{
			}
		case 25:
			if(tailed == ONE_TAILED){
				if (t > 1.708141){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.059539 || t > 2.059539){
					return false;
				}
			}else{
			}
		case 26:
			if(tailed == ONE_TAILED){
				if (t > 1.705618){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.055529 || t > 2.055529){
					return false;
				}
			}else{
			}
		case 27:
			if(tailed == ONE_TAILED){
				if (t > 1.703288){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.051831 || t > 2.051831){
					return false;
				}
			}else{
			}
		case 28:
			if(tailed == ONE_TAILED){
				if (t > 1.701131){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.048407 || t > 2.048407){
					return false;
				}
			}else{
			}
		case 29:
			if(tailed == ONE_TAILED){
				if (t > 1.699127){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.04523 || t > 2.04523){
					return false;
				}
			}else{
			}
		case 30:
			if(tailed == ONE_TAILED){
				if (t > 1.697261){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.042272 || t > 2.042272){
					return false;
				}
			}else{
			}
		}
	}else if (alpha == 0.025){
		switch(df){
		case 1:
			if(tailed == ONE_TAILED){
				if (t > 12.7062){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -25.4517 || t > 25.4517){
					return false;
				}
			}else{
			}
		case 2:
			if(tailed == ONE_TAILED){
				if (t > 4.302653){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -6.205347 || t > 6.205347){
					return false;
				}
			}else{
			}
		case 3:
			if(tailed == ONE_TAILED){
				if (t > 3.182446){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -4.176535 || t > 4.176535){
					return false;
				}
			}else{
			}
		case 4:
			if(tailed == ONE_TAILED){
				if (t > 2.776298){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -3.494943 || t > 3.494943){
					return false;
				}
			}else{
			}
		case 5:
			if(tailed == ONE_TAILED){
				if (t > 2.570543){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -3.163262 || t > 3.163262){
					return false;
				}
			}else{
			}
		case 6:
			if(tailed == ONE_TAILED){
				if (t > 2.446899){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.968647 || t > 2.968647){
					return false;
				}
			}else{
			}
		case 7:
			if(tailed == ONE_TAILED){
				if (t > 2.364619){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.841229 || t > 2.841229){
					return false;
				}
			}else{
			}
		case 8:
			if(tailed == ONE_TAILED){
				if (t > 2.306002){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.751517 || t > 2.751517){
					return false;
				}
			}else{
			}
		case 9:
			if(tailed == ONE_TAILED){
				if (t > 2.262156){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.685007 || t > 2.685007){
					return false;
				}
			}else{
			}
		case 10:
			if(tailed == ONE_TAILED){
				if (t > 2.228138){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.633765 || t > 2.633765){
					return false;
				}
			}else{
			}
		case 11:
			if(tailed == ONE_TAILED){
				if (t > 2.200985){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.593092 || t > 2.593092){
					return false;
				}
			}else{
			}
		case 12:
			if(tailed == ONE_TAILED){
				if (t > 2.178813){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.560032 || t > 2.560032){
					return false;
				}
			}else{
			}
		case 13:
			if(tailed == ONE_TAILED){
				if (t > 2.160369){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.532637 || t > 2.532637){
					return false;
				}
			}else{
			}
		case 14:
			if(tailed == ONE_TAILED){
				if (t > 2.144787){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.509569 || t > 2.509569){
					return false;
				}
			}else{
			}
		case 15:
			if(tailed == ONE_TAILED){
				if (t > 2.131449){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.48988 || t > 2.48988){
					return false;
				}
			}else{
			}
		case 16:
			if(tailed == ONE_TAILED){
				if (t > 2.119905){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.472878 || t > 2.472878){
					return false;
				}
			}else{
			}
		case 17:
			if(tailed == ONE_TAILED){
				if (t > 2.109816){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.458051 || t > 2.458051){
					return false;
				}
			}else{
			}
		case 18:
			if(tailed == ONE_TAILED){
				if (t > 2.100922){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.445006 || t > 2.445006){
					return false;
				}
			}else{
			}
		case 19:
			if(tailed == ONE_TAILED){
				if (t > 2.093024){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.43344 || t > 2.43344){
					return false;
				}
			}else{
			}
		case 20:
			if(tailed == ONE_TAILED){
				if (t > 2.085963){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.423116 || t > 2.423116){
					return false;
				}
			}else{
			}
		case 21:
			if(tailed == ONE_TAILED){
				if (t > 2.079614){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.413845 || t > 2.413845){
					return false;
				}
			}else{
			}
		case 22:
			if(tailed == ONE_TAILED){
				if (t > 2.073873){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.405473 || t > 2.073873){
					return false;
				}
			}else{
			}
		case 23:
			if(tailed == ONE_TAILED){
				if (t > 2.068658){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.397875 || t > 2.397875){
					return false;
				}
			}else{
			}
		case 24:
			if(tailed == ONE_TAILED){
				if (t > 2.063899){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.390949 || t > 2.390949){
					return false;
				}
			}else{
			}
		case 25:
			if(tailed == ONE_TAILED){
				if (t > 2.059539){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.38461 || t > 2.38461){
					return false;
				}
			}else{
			}
		case 26:
			if(tailed == ONE_TAILED){
				if (t > 2.055529){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.378786 || t > 2.378786){
					return false;
				}
			}else{
			}
		case 27:
			if(tailed == ONE_TAILED){
				if (t > 2.051831){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.373417 || t > 2.373417){
					return false;
				}
			}else{
			}
		case 28:
			if(tailed == ONE_TAILED){
				if (t > 2.048407){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.368452 || t > 2.368452){
					return false;
				}
			}else{
			}
		case 29:
			if(tailed == ONE_TAILED){
				if (t > 2.04523){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.363846 || t > 2.363846){
					return false;
				}
			}else{
			}
		case 30:
			if(tailed == ONE_TAILED){
				if (t > 2.042272){
					return false;
				}
			}else if(tailed == TWO_TAILED){
				if (t < -2.359562 || t > 2.359562){
					return false;
				}
			}else{
			}
		}

	}else{
	}

	return true;
}

template <class T, int size> double MyStatistics<T, size>::tWhereToStartCounting(int df){
	switch(df){
	case 1:
		return 564.188700;
	case 2:
		return 99.990000;
	case 3:
		return 42.612100;
	case 4:
		return 25.974900;
	case 5:
		return 18.895400;
	case 6:
		return 15.172100;
	case 7:
		return 12.937000;
	case 8:
		return 11.469600;
	case 9:
		return 10.442100;
	case 10:
		return 9.687300;
	case 11:
		return 9.111800;
	case 12:
		return 8.659700;
	case 13:
		return 8.296100;
	case 14:
		return 7.997700;
	case 15:
		return 7.748700;
	case 16:
		return 7.538000;
	case 17:
		return 7.357500;
	case 18:
		return 7.201200;
	case 19:
		return 7.064700;
	case 20:
		return 6.944400;
	case 21:
		return 6.837600;
	case 22:
		return 6.742300;
	case 23:
		return 6.656600;
	case 24:
		return 6.579200;
	case 25:
		return 6.509000;
	case 26:
		return 6.445000;
	case 27:
		return 6.386500;
	case 28:
		return 6.332700;
	case 29:
		return 6.283100;
	case 30:
		return 6.237300;
	}
	return 0;
}

/*
ciSquared stuff*******
*/

// colExpeced & colActual are indexes

template <class T, int size> bool MyStatistics<T, size>::chiSquaredTest(int colExpeced, int colActual, double alpha){
	double df = MyData<T, size>::recordCount() - 1;
	double c = chiSquared(colExpeced, colActual);
	double area = chiSquaredIntegral(c, df);
	if (area >= alpha){
		return true;//fail to reject null hypothesis
	}else{
		return false;
	}
}


template <class T, int size> double MyStatistics<T, size>::chiSquared(int colExpeced, int colActual){
	typename MyData<T, size>::iterator i = MyData<T, size>::begin();
	double p_value = 0;
	while (i != MyData<T, size>::end()){
		myArray temp = *i;
		p_value += (double)pow(temp[colActual] - temp[colExpeced], 2)/temp[colExpeced];
		++i;
	}

	return p_value;
}

template <class T, int size> double MyStatistics<T, size>::chiSquaredProbabilityDensity(double df, double x){
	const double E = 2.71828182845904523536;//there's probably a better place for this.
	double tg = tgamma(df/2);
	double tf = (double)1/(pow(2,df/2)*tg);
	double ts = pow(x,(df/2)-1);
	double tt = pow(E,-x/2);
	double t = tf*ts*tt;
	return t;
}


template <class T, int size> double MyStatistics<T, size>::chiSquaredIntegral(double cX, double df, double x){
	double area = 0;
	double dx = (double)1/pow(10, 4);
	double n = 30.0;
	while (n >= cX){
		area += Ctrapezoid( n, dx, df);
		n -= dx;
	}
	return area;
}

template <class T, int size> double MyStatistics<T, size>::Ctrapezoid(double t, double dx, int df){
	double a = chiSquaredProbabilityDensity(df, t);
	double b = chiSquaredProbabilityDensity(df, t - dx);
	return ((double)(a + b)/2.0 * dx);
}






template <class T, int size> int MyStatistics<T, size>::exceptionHandler(int error){
	switch(error){
	case ERROR_OUT_OF_RANGE:
		cerr << "Error: " << to_string(error) << " OUT OF BOUNDS" << '\n';
		break;
	}

	return error;
}

/*
This stuff was useful when I was trying to figure out how the t score stuff is calculated.
I don't know if it's still useful.
*/

/*
* This is where the switch statement in halfAreaUnderTDistribution() came from.
*/
template <class T, int size> void MyStatistics<T, size>::areaForAllCurves(){
	cerr << "\tswitch(df){" << '\n';
	for (int x = 1; x <= 30; x++){
		double area = halfTotalAreaUnderCurve(x);
		cerr << "\tcase " << to_string(x) << ":" << '\n';
		cerr << "\t\treturn " << to_string(area) << ";" << '\n';
	}
	cerr << "\t}" << '\n';
}

/*
* This will return half the area under the curve for whatever degrees of freedom you give it.
*/
template <class T, int size> double MyStatistics<T, size>::halfTotalAreaUnderCurve(int df){
	double area = 0;
	double dx = (double)1/pow(10, 4);
	double n = 0.0;

	while( ProbabilityDistributionFormula(n, df) >= 0.0000001){
		area += Ttrapezoid( n, dx);
		n += dx;
	}

	return area;
}



/*
*	multiply by 2 to get full area under T distribution curve for given degrees of freedom.
*/
template <class T, int size> double MyStatistics<T, size>::halfAreaUnderTDistribution(int df){
	switch(df){
	case 1:
		return 0.499436;
	case 2:
		return 0.496817;
	case 3:
		return 0.492531;
	case 4:
		return 0.487752;
	case 5:
		return 0.483170;
	case 6:
		return 0.479050;
	case 7:
		return 0.475444;
	case 8:
		return 0.472318;
	case 9:
		return 0.469609;
	case 10:
		return 0.467257;
	case 11:
		return 0.465205;
	case 12:
		return 0.463404;
	case 13:
		return 0.461816;
	case 14:
		return 0.460405;
	case 15:
		return 0.459147;
	case 16:
		return 0.458018;
	case 17:
		return 0.457000;
	case 18:
		return 0.456079;
	case 19:
		return 0.455241;
	case 20:
		return 0.454476;
	case 21:
		return 0.453775;
	case 22:
		return 0.453131;
	case 23:
		return 0.452536;
	case 24:
		return 0.451986;
	case 25:
		return 0.451476;
	case 26:
		return 0.451002;
	case 27:
		return 0.450560;
	case 28:
		return 0.450147;
	case 29:
		return 0.449760;
	case 30:
		return 0.449397;
	}

	return -1.0;//create an error
}



