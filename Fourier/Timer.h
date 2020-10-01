//Timer.h : main header file

#if !defined(_TIMER_H_)
#define _TIMER_H_

#include <winbase.h>

class CTimer
{
private:
	LARGE_INTEGER litmp;
	LONGLONG llQFreq,llQPart1,llQPart2;
	double fFreq;

public:
	CTimer()
	{
		QueryPerformanceFrequency(&litmp);
		llQFreq=litmp.QuadPart;
		fFreq=(double)llQFreq;
	}

	~CTimer() {;}

public:
	inline void Start()
	{
		QueryPerformanceCounter(&litmp);
		llQPart1=litmp.QuadPart;
	}

	inline double End()
	{
		QueryPerformanceCounter(&litmp);
		llQPart2=litmp.QuadPart;
		double ftime=((double)(llQPart2-llQPart1))/fFreq;
		llQPart1=llQPart2;
		return ftime;
	}
};

#endif //!defined(_TIMER_H_)
