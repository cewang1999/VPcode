/*Field is a container for mesh node data division by volume*/
#ifndef _FIELD_H
#define _FIELD_H
#include <vector>
#include <ostream>

// *** vec2 ****
template <typename T>
struct vec2 {
	vec2 (const T u, const T v) : d{u,v} {}
	vec2 (const T a[2]) : d{a[0],a[1]} {}
	vec2 (): d{0,0} {}
	T& operator[](int i) {return d[i];}
	T operator()(int i) const {return d[i];}
	vec2<T>& operator=(double s) {d[0]=s;d[1]=s;return (*this);}
	vec2<T>& operator+=(vec2<T> o) {d[0]+=o[0];d[1]+=o[1];return(*this);}
	vec2<T>& operator-=(vec2<T> o) {d[0]-=o[0];d[1]-=o[1];return(*this);}
	vec2<T> operator/(double s) {vec2<T>o; o[0]=d[0]/s;o[1]=d[1]/s;return o;}
    vec2<T> operator/=(double s) {d[0]/=s;d[1]/=s;return (*this);}

	//dot product of two vectors
	friend T dot(const vec2<T> &v1, const vec2<T> &v2) {
		T s=0;	for (int i=0;i<2;i++) s+=v1(i)*v2(i);
		return s;	}

	//vector magnitude
	friend T mag(const vec2<T> &v) {return sqrt(dot(v,v));}

	//unit vector
	friend vec2<T> unit(const vec2<T> &v) {return vec2(v)/mag(v);}

protected:
	T d[2];
};

//vec2-vec2 operations
template<typename T>	//addition of two vec3s
vec2<T> operator+(const vec2<T>& a, const vec2<T>& b) {
	return vec2<T> (a(0)+b(0),a(1)+b(1));	}
template<typename T>	//subtraction of two vec2s
vec2<T> operator-(const vec2<T>& a, const vec2<T>& b) {
	return vec2<T> (a(0)-b(0),a(1)-b(1));	}
template<typename T>	//element-wise multiplication of two vec2s
vec2<T> operator*(const vec2<T>& a, const vec2<T>& b) {
	return vec2<T> (a(0)*b(0),a(1)*b(1));	}
template<typename T>	//element wise division of two vec3s
vec2<T> operator/(const vec2<T>& a, const vec2<T>& b) {
	return vec2<T> (a(0)/b(0),a(1)/b(1));	}

//vec2 - scalar operations
template<typename T>		//scalar multiplication
vec2<T> operator*(const vec2<T> &a, T s) {
	return vec2<T>(a(0)*s, a(1)*s);}
template<typename T>		//scalar multiplication 2
vec2<T> operator*(T s,const vec2<T> &a) {
	return vec2<T>(a(0)*s, a(1)*s);}

//output
template<typename T>	//ostream output
std::ostream& operator<<(std::ostream &out, vec2<T>& v) {
	out<<v[0]<<" "<<v[1]<<" 0";	//paraview does not support 2-component arrays
	return out;
}

using double2 = vec2<double>;
using int2 = vec2<int>;


// *** Field_ ************
template <typename T>
class Field_
{
public:
	
	/*constructor*/
	Field_(int ni, int nj) :
	ni{ni}, nj{nj}
	{
		//allocate memory for a 2D array
		data = new T*[ni];
		for (int i=0;i<ni;i++)
			data[i] = new T[nj];

		clear();
	}

	//another constructor taking an int2
	Field_(int2 nn) : Field_(nn[0],nn[1]) {};

	//copy constructor
	Field_(const Field_ &other):
	Field_{other.ni,other.nj} {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				data[i][j] = other(i,j);
	}
	
	//move constructor
	Field_(Field_ &&other):	ni{other.ni},nj{other.nj} {
			data = other.data;	//steal the data
			other.data = nullptr;	//invalidate
	}

	//move assignment operator
	Field_& operator = (Field_ &&f) {data=f.data;
				f.data=nullptr; return *this;}

	//destructor: release memory
	~Field_() {
		//don't do anything if data is not allocated (or was moved away)
		if (data==nullptr) return;
		for (int i=0;i<ni;i++)	delete[] data[i];
		delete[] data;
	}

	//overloaded operator [] to allow direct access to data
	T* operator[] (int i) {return data[i];}

	/*returns data[i][j] marked as const to signal no data change*/
	T operator() (int i, int j) const {return data[i][j];}

	/*sets all values to some scalar*/
	void operator =(double s) {
		for (int i=0;i<ni;i++)
		  for (int j=0;j<nj;j++)
		   	data[i][j] = s;
	  }

	/*performs element by element division by another field*/
	void operator /= (const Field_ &other) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
			{
				if (other.data[i][j]!=0)
				  data[i][j] /= other.data[i][j];
				else
				  data[i][j] = 0;
			}  
	}

	/*increments values by data from another field*/
	Field_& operator += (const Field_ &other) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				data[i][j]+=other(i,j);
		return (*this);
	}

	/*performs element by element division by another field*/
	Field_& operator *= (double s) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				data[i][j]*=s;
		return (*this);
	}

	//multiplication operator, returns f*s
	friend Field_<T> operator*(double s, const Field_<T>&f) {
		Field_<T> r(f);
		return r*=s;
	}

	//division of a field by a field of doubles
	friend Field_<T> operator*(const Field_<T>&f1, const Field_<T>&f2) {
		Field_<T> r(f1);
		for (int i=0;i<f1.ni;i++)
			for (int j=0;j<f1.nj;j++)
				r[i][j] = f1(i,j)*f2(i,j);
		return r;
	}

	//division of a field by a field of doubles
	friend Field_<T> operator/(const Field_<T>&f, const Field_<double>&d) {
		Field_<T> r(f);
		for (int i=0;i<f.ni;i++)
			for (int j=0;j<f.nj;j++)
			{
				if (d(i,j)!=0)	//check for div by zero
					r[i][j] = f(i,j)/d(i,j);
				else
					r[i][j] = 0;
			}
		return r;
	}


	/*sets all data to zero*/
	void clear() {(*this)=0;}


	const int ni,nj;	//allocated dimensions

protected:
	T **data;	/*data held by this field*/
};

/*writes out data to a file stream*/
template<typename T>
std::ostream& operator<<(std::ostream &out, Field_<T> &f)
{
	for (int j=0;j<f.nj;j++)
		for (int i=0;i<f.ni;i++) out<<f.data[i][j]<<" ";
	return out;
}

//some typedefs
using Field  = Field_<double>;
using FieldI = Field_<int>;
using dvector = std::vector<double>;
using dvectorI = std::vector<int>;
using namespace std;
#endif
