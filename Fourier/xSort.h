//xSort.h : main header file

#if !defined(_XSORT_H_)
#define _XSORT_H_

template<class T> class CxSort//≈≈–ÚÀ„∑®: BubbleSort,ShellSort,MergeSort
{
public:
	CxSort() {;}
	~CxSort() {;}

private:
	inline void Merge(T* source,T* result,int size,int n);

public:
	inline void BubbleSort(int nLength,T* tData);
	inline void ShellSort(int nLength,T* ptData);
	inline bool MergeSort(int nLength,T* ptData);
};

template<class T> void CxSort<T>::BubbleSort(int nLength,T* tData)
{   
	for(int i=0;i<nLength;i++)//Make a pass through the array for each element
	{
		for(int j=1;j<nLength-i;j++)//Go through the array beginning to end
		{
			if(tData[j-1]>tData[j])//If the the first number is greater,swap it 
			{
				T x=tData[j]; tData[j]=tData[j-1]; tData[j-1]=x;
			}
		}
	}
}

template<class T> void CxSort<T>::ShellSort(int nLength,T* tData)
{
	int i,j,increment=3;
	T temp;
	while(increment>0)
	{
		for(i=0;i<nLength;i++)
		{
			j=i; temp=tData[i];
			while((j>=increment) && (tData[j-increment]>temp))
			{
				tData[j]=tData[j-increment];
				j=j-increment;
			}
			tData[j]=temp;
		}
		if(increment/2!=0) increment=increment/2;
		else if(1==increment) increment=0;
		else increment=1;
	}
}

template<class T> void CxSort<T>::Merge(T* source,T* result,int size,int n)
{
	int lb1=0,lb2,ub1,ub2,p=0,i,j;
	while((lb1+size)<n)
	{
		lb2=lb1+size; ub1=lb2-1;
		if((lb2+size-1)>n) ub2=n-1;
		else ub2=lb2+size-1;
		i=lb1; j=lb2;
		while((i<=ub1)&&(j<=ub2))
		{
			if(source[i]<=source[j]) result[p++]=source[i++];
			else result[p++]=source[j++];
		}
		while(i<=ub1) result[p++]=source[i++];
		while(j<=ub2) result[p++]=source[j++];
		lb1=ub2+1;
	}
	i=lb1;
	while(p<n) result[p++]=source[i++];
}

template<class T> bool CxSort<T>::MergeSort(int nLength,T* ptData)
{
	int s=1,n=nLength;
	T* temp=new T[n];

	if(temp==NULL) return false;

	while(s<n)
	{
		Merge(ptData,temp,s,n); s*=2;
		Merge(temp,ptData,s,n); s*=2;
	}

	delete[] temp;
	return true;
}

#endif //!defined(_XSORT_H_)
