#ifndef _DATA_H_
#define _DATA_H_



class Data
{
  public:
  	
  	Data();
  	virtual ~Data();
  	
	void setinputNumber(int iN);
	int getinputNumber();
    
    void setDataPointer(double* d);
    double* getDataPointer();
    
    void setFZCount(int fzcount);
    int getFZCount();
    
    void setPointGroup(int pgnum);
    int getPointGroup();
    
    void setNumSteps(int nsteps);
    int getNumSteps();
    
    void execute();

  private:
    int m_PointGroup;
    int m_NumSteps; 
    int m_FZCount;
    double* m_Data;
};

#endif /* _DATA_H_ */