#include <Solver/CppDASSL/CppDASSL.h>
#include <Core/Math/Functions.h>
//#include <Core/Math/ILapack.h>
#include <fstream>
#include <iomanip>

#if defined(USE_OPENMP)
#include "omp.h"
#include <Core/Utils/numeric/bindings/umfpack/umfpack.hpp>
#include <Core/Utils/numeric/bindings/ublas/vector.hpp>
#include <Core/Utils/numeric/bindings/ublas.hpp>
#include <boost/numeric/ublas/io.hpp>

CppDASSL::CppDASSL(IMixedSystem* system, ISolverSettings* settings)
    : SolverDefaultImplementation(system, settings),
      _cppdasslsettings(dynamic_cast<ISolverSettings*>(_settings)),
      _continuous_system(),
      _time_system()
/*      _cvodeMem(NULL),
      _z(NULL),
      _zInit(NULL),
      _zWrite(NULL),
      _dimSys(0),
      _cv_rt(0),
      _outStps(0),
      _locStps(0),
      _idid(0),
      _hOut(0.0),
      _tOut(0.0),
      _tZero(0.0),
      _zeroSign(NULL),
      _absTol(NULL),
      _cvode_initialized(false),
      _tLastEvent(0.0),
      _event_n(0),
      _properties(NULL),
      _continuous_system(NULL),
      _event_system(NULL),
      _mixed_system(NULL),
      _time_system(NULL),
    _delta(NULL),
    _ysave(NULL) */
{
  _data = ((void*) this);
}

CppDASSL::~CppDASSL()
{
  if (_y)
    delete [] _y;
  if (_yp)
    delete [] _yp;
  if(_continuous_system) {
    delete [] _continuous_system;
  }
}

void CppDASSL::initialize()
{
    SolverDefaultImplementation::initialize();
    IContinuous *continuous_system = dynamic_cast<IContinuous*>(_system);
    ITime *time_system =  dynamic_cast<ITime*>(_system);
    int numThreads=_cppdasslsettings->getGlobalSettings()->getSolverThreads();
    dasslSolver.setNumThreads(numThreads);
    _continuous_system = new IContinuous*[numThreads];
    _time_system = new ITime*[numThreads];
    _continuous_system[0] = continuous_system;
    _time_system[0] = time_system;
    _event_system = dynamic_cast<IEvent*>(_system);
    _mixed_system = dynamic_cast<IMixedSystem*>(_system);
    for(int i = 1; i < numThreads; i++)
    {
        IMixedSystem* clonedSystem = _system->clone();
        _continuous_system[i] = dynamic_cast<IContinuous*>(clonedSystem);
        _time_system[i] = dynamic_cast<ITime*>(clonedSystem);
        dynamic_cast<ISystemInitialization*>(clonedSystem)->initialize();
    }


    _dimSys = _continuous_system[0]->getDimContinuousStates();
    _dimZeroFunc = _event_system->getDimZeroFunc();
    _y = new double[_dimSys];
    _yp = new double[_dimSys];


    _continuous_system[0]->evaluateAll(IContinuous::ALL);
    _continuous_system[0]->getContinuousStates(_y);
// begin analyzation mode
    _continuous_system[0]->setContinuousStates(_y);
    _continuous_system[0]->evaluateODE(IContinuous::ALL);    // vxworksupdate
    _continuous_system[0]->getRHS(_yp);

    int _countnz=0;
    double delta=1e-6;
    double* _yphelp=new double[_dimSys];
    _time_system[0]->setTime(_tCurrent);
    for(int i=0; i<_dimSys; ++i) {
        double _ysave;
        _ysave=_y[i];
        _y[i]+=delta;
        _continuous_system[0]->setContinuousStates(_y);
        _continuous_system[0]->evaluateODE(IContinuous::ALL);    // vxworksupdate
        _continuous_system[0]->getRHS(_yphelp);
        for(int j=0; j<_dimSys; ++j) {
            if(_yphelp[j]-_yp[j]>1e-12) {
                _countnz++;
            }
        }
        if(_countnz>_dimSys*_dimSys/100) {
            dasslSolver.setSparse(false);
            break;
        }
        _y[i]=_ysave;
    }
    if(_countnz<_dimSys*_dimSys/100) {
        dasslSolver.setSparse(true);
    }
    dasslSolver.setDenseOutput(true);
    delete [] _yphelp;
}



void CppDASSL::solve(const SOLVERCALL action)
{
    double t;
    if ((action & RECORDCALL) && (action & FIRST_CALL)) {
        initialize();
        return;
    }

    if (action & RECALL)
    {
        writeToFile(0, _tCurrent, _h);
        _continuous_system[0]->getContinuousStates(_y);
    }
  // Initialization phase
    t=_tCurrent;
    _time_system[0]->setTime(t);
    _continuous_system[0]->setContinuousStates(_y);
    _continuous_system[0]->evaluateODE(IContinuous::ALL);    // vxworksupdate
    _continuous_system[0]->getRHS(_yp);
    SolverDefaultImplementation::writeToFile(0, t, _h);
    bool state_selection = stateSelection();
    if(!_dimZeroFunc) {
        int idid=dasslSolver.solve(&res,_dimSys,t,&_y[0],&_yp[0],_tEnd,_data,NULL,NULL,NULL,_dimZeroFunc,NULL,false);
        _continuous_system[0]->stepCompleted(t);
        if (state_selection) {
            _continuous_system[0]->getContinuousStates(_y);
            _continuous_system[0]->setContinuousStates(_y);
            _continuous_system[0]->evaluateODE(IContinuous::ALL);    // vxworksupdate
            _continuous_system[0]->getRHS(_yp);
        }
        while(idid==-1 || idid==1) {
            _time_system[0]->setTime(t);
            _continuous_system[0]->setContinuousStates(_y);
            _continuous_system[0]->evaluateAll(IContinuous::ALL);
            SolverDefaultImplementation::writeToFile(0, t, _h);
            idid=dasslSolver.solve(&res,_dimSys,t,&_y[0],&_yp[0],_tEnd,_data,NULL,NULL,NULL,_dimZeroFunc,NULL,true);
            _continuous_system[0]->stepCompleted(t);
            if (state_selection) {
                _continuous_system[0]->getContinuousStates(_y);
                _continuous_system[0]->setContinuousStates(_y);
                _continuous_system[0]->evaluateODE(IContinuous::ALL);    // vxworksupdate
                _continuous_system[0]->getRHS(_yp);
            }
        }
    } else {
        if(_jroot) delete [] _jroot;
        _jroot=new int[_dimZeroFunc];
        int idid=dasslSolver.solve(&res,_dimSys,t,&_y[0],&_yp[0],_tEnd,_data,NULL,NULL,&zeroes,_dimZeroFunc,_jroot,false);
        _continuous_system[0]->stepCompleted(t);
        if (state_selection) {
            _continuous_system[0]->getContinuousStates(_y);
            _continuous_system[0]->setContinuousStates(_y);
            _continuous_system[0]->evaluateODE(IContinuous::ALL);    // vxworksupdate
            _continuous_system[0]->getRHS(_yp);
        }
        if(idid==5) {
            _time_system[0]->setTime(t);
            _continuous_system[0]->setContinuousStates(_y);
            _continuous_system[0]->evaluateAll(IContinuous::ALL);
            SolverDefaultImplementation::writeToFile(0, t, _h);
            for (int i = 0; i < _dimZeroFunc; i++) _events[i] = bool(_jroot[i]);

            if (_mixed_system->handleSystemEvents(_events))
              {
                _continuous_system[0]->getContinuousStates(_y);
                _continuous_system[0]->setContinuousStates(_y);
                _continuous_system[0]->evaluateODE(IContinuous::ALL);    // vxworksupdate
                _continuous_system[0]->getRHS(_yp);
              }

        }
        while(idid==-1 || idid==5 || idid==1) {
            _time_system[0]->setTime(t);
            _continuous_system[0]->setContinuousStates(_y);
            _continuous_system[0]->evaluateAll(IContinuous::ALL);
            SolverDefaultImplementation::writeToFile(0, t, _h);
            if(idid==5) {
                idid=dasslSolver.solve(&res,_dimSys,t,&_y[0],&_yp[0],_tEnd,_data,NULL,NULL,&zeroes,_dimZeroFunc,_jroot,false);
            } else {
                idid=dasslSolver.solve(&res,_dimSys,t,&_y[0],&_yp[0],_tEnd,_data,NULL,NULL,&zeroes,_dimZeroFunc,_jroot,true);
            }
            if(_continuous_system[0]->stepCompleted(t)) break;
            if (state_selection) {
                _continuous_system[0]->getContinuousStates(_y);
                _continuous_system[0]->setContinuousStates(_y);
                _continuous_system[0]->evaluateODE(IContinuous::ALL);    // vxworksupdate
                _continuous_system[0]->getRHS(_yp);
            }
;
            if(idid==5) {
                _time_system[0]->setTime(t);
                _continuous_system[0]->setContinuousStates(_y);
                _continuous_system[0]->evaluateAll(IContinuous::ALL);
                SolverDefaultImplementation::writeToFile(0, t, _h);
                for (int i = 0; i < _dimZeroFunc; i++) _events[i] = bool(_jroot[i]);

                if (_mixed_system->handleSystemEvents(_events))
                  {
                    _continuous_system[0]->getContinuousStates(_y);
                    _continuous_system[0]->setContinuousStates(_y);
                    _continuous_system[0]->evaluateODE(IContinuous::ALL);    // vxworksupdate
                    _continuous_system[0]->getRHS(_yp);
                  }

            }
        }

    }
    _tCurrent=_tEnd;
    _time_system[0]->setTime(_tCurrent);
    _continuous_system[0]->setContinuousStates(_y);
    _continuous_system[0]->evaluateAll(IContinuous::ALL);
    SolverDefaultImplementation::writeToFile(0, t, _h);
    _solverStatus = ISolver::DONE;
}



void CppDASSL::writeCppDASSLOutput(const double &time, const double &h, const int &stp)
{
  /*
  if (stp > 0)
  {
    if (_cvodesettings->getDenseOutput())
    {
      _bWritten = false;
      double *oldValues = NULL;

      //We have to find all output-points within the last solver step
      while (_tLastWrite + dynamic_cast<ISolverSettings*>(_cvodesettings)->getGlobalSettings()->gethOutput() <= time)
      {
        if (!_bWritten)
        {
          //Rescue the calculated derivatives
          oldValues = new double[_continuous_system->getDimRHS()];
          _continuous_system->getRHS(oldValues);
        }
        _bWritten = true;
        _tLastWrite = _tLastWrite + dynamic_cast<ISolverSettings*>(_cvodesettings)->getGlobalSettings()->gethOutput();
        //Get the state vars at the output-point (interpolated)
        _idid = CVodeGetDky(_cvodeMem, _tLastWrite, 0, _CV_yWrite);
        _time_system->setTime(_tLastWrite);
        _continuous_system->setContinuousStates(NV_DATA_S(_CV_yWrite));
        _continuous_system->evaluateAll(IContinuous::CONTINUOUS);
        SolverDefaultImplementation::writeToFile(stp, _tLastWrite, h);
      }      //end if time -_tLastWritten
      if (_bWritten)
      {
        _time_system->setTime(time);
        _continuous_system->setContinuousStates(_z);
        _continuous_system->setRHS(oldValues);
        delete[] oldValues;
        //_continuous_system->evaluateAll(IContinuous::CONTINUOUS);
      }
      else if (time == _tEnd && _tLastWrite != time)
      {
        _idid = CVodeGetDky(_cvodeMem, time, 0, _CV_y);
        _time_system->setTime(time);
        _continuous_system->setContinuousStates(NV_DATA_S(_CV_y));
        _continuous_system->evaluateAll(IContinuous::CONTINUOUS);
        SolverDefaultImplementation::writeToFile(stp, _tEnd, h);
      }
    }
    else
      SolverDefaultImplementation::writeToFile(stp, time, h);
  }
  */
}

int CppDASSL::res(const double* t, const double* y, const double* yprime, double* cj, double* delta, int* ires, void *par) {
    return ((CppDASSL*) par)->calcFunction(t, y, yprime, cj, delta, ires);
}

int CppDASSL::calcFunction(const double* t, const double* y, const double* yprime, double* cj, double* delta, int* ires) {
    int numThread=omp_get_thread_num();
    _time_system[numThread]->setTime(*t);
    _continuous_system[numThread]->setContinuousStates(y);
    _continuous_system[numThread]->evaluateODE(IContinuous::ALL);    // vxworksupdate
    _continuous_system[numThread]->getRHS(delta);
    for(int i=0; i<_dimSys; ++i) delta[i]=yprime[i]-delta[i];
    return 0;
}

int CppDASSL::zeroes(const int* NEQ, const double* T, const double* Y, const double* YP, int* NRT, double* RVAL, void *par)
{
  ((CppDASSL*) par)->giveZeroVal(*T, Y, RVAL);

  return (0);
}


void CppDASSL::giveZeroVal(const double &t, const double *y, double *zeroValue)
{
  _time_system[0]->setTime(t);
  _continuous_system[0]->setContinuousStates(y);

  // System aktualisieren
  _continuous_system[0]->evaluateZeroFuncs(IContinuous::DISCRETE);
  _event_system->getZeroFunc(zeroValue);
}

bool CppDASSL::stateSelection()
{
  return SolverDefaultImplementation::stateSelection();
}

void CppDASSL::setTimeOut(unsigned int time_out)
  {
       SimulationMonitor::setTimeOut(time_out);
  }
 void CppDASSL::stop()
  {
       SimulationMonitor::stop();
  }


//void CppDASSL::setcycletime(double cycletime){}
void CppDASSL::writeSimulationInfo(){}
int CppDASSL::reportErrorMessage(std::ostream& messageStream) {
    return 0;
}
#endif
