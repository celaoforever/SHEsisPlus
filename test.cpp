#include <iostream>
#include <stdio.h>
template  <class T>
class array1D {
public:
	array1D(T a[],int msize):a(a),msize(msize),isStatic(true){};
	~array1D(){
		if(!isStatic){
			delete[] a;
			a=0;
		}
	}
	int size(){return msize;};
	void SetIsStatic(bool b)
	{this->isStatic=b;};
	T& operator [] (const int idx){
		if(idx<msize)
			return a[idx];
		else
			std::cout<<"error";
	}
private:
	T* a;
	int msize;
	bool isStatic;
};

int main(){
int array[]={1,1,1};
array1D<int> a({1,1,1},3);
for(int i=0;i<a.size();i++)
{
	std::cout<<a[i]<<",";
};
return 0;
};
