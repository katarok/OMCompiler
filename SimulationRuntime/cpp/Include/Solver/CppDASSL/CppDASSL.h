
#pragma once

#include "FactoryExport.h"
#include <Core/ModelicaDefine.h>
#include <Core/Modelica.h>
#include <Core/Solver/SolverDefaultImplementation.h>
#include <Core/Utils/extension/measure_time.hpp>

int ddaskr_ (
    int (*RES)(const double* t,const double* y,const double* yprime, double* cj, double* delta, int* ires, void *),
    int* NEQ,
    double* T,
    double* Y,
    double* YPRIME,
    double* TOUT,
    int* INFO,
    double* RTOL,
    double* ATOL,
    int* IDID,
    double* RWORK,
    int* LRW,
    int* IWORK,
    int* LIW,
    void* PAR,
    int (*JAC)(double *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, void *),
    int (*PSOL)(int *NEQ, double* T, double* Y, double* YPRIME, double* SAVR, double* WK, double* CJ, double* WGHT, double * WP, int* IWP, double* B, double* EPLIN, int* IER, void *),
    int (*RT)(int* NEQ,double* T,double* Y,double* YP,int* NRT, double* RVAL, void *),
    int* NRT,
    int* JROOT
);

#if defined(USE_OPENMP)
/*****************************************************************************/
// Peer
// BDF-Verfahren für steife und nicht-steife ODEs
// Dokumentation siehe offizielle Peer Doku

/*****************************************************************************
Copyright (c) 2014, IWR TU Dresden, All rights reserved
*****************************************************************************/

class CppDASSL
  : public ISolver,  public SolverDefaultImplementation
{
public:

  CppDASSL(IMixedSystem* system, ISolverSettings* settings);

  virtual ~CppDASSL();

  // geerbt von Object (in SolverDefaultImplementation)
  //---------------------------------------
  /// Spezielle Solvereinstellungen setzten (default oder user defined)
  virtual void initialize();


  // geerbt von ISolver
  //---------------------------------------
  /// Setzen der Startzeit für die numerische Lösung
  virtual void setStartTime(const double& time)
  {
    SolverDefaultImplementation::setStartTime(time);
  };

  /// Setzen der Endzeit für die numerische Lösung
  virtual void setEndTime(const double& time)
  {
    SolverDefaultImplementation::setEndTime(time);
  };

  /// Setzen der initialen Schrittweite (z.B. auch nach Nullstelle)
  virtual void setInitStepSize(const double& stepSize)
  {
    SolverDefaultImplementation::setInitStepSize(stepSize);
  };

  /// Berechung der numerischen Lösung innerhalb eines gegebenen Zeitintervalls
  virtual void solve(const SOLVERCALL command = UNDEF_CALL);

  /// Liefert den Status des Solvers nach Beendigung der Simulation
  virtual ISolver::SOLVERSTATUS getSolverStatus()
  {
    return (SolverDefaultImplementation::getSolverStatus());
  };

  //// Ausgabe von statistischen Informationen (wird vom SimManager nach Abschluß der Simulation aufgerufen)
  virtual void writeSimulationInfo();


  virtual int reportErrorMessage(std::ostream& messageStream);
  virtual bool stateSelection();
  virtual void setTimeOut(unsigned int time_out);
    virtual void stop();
private:





  // Nulltellenfunktion
  void writeCppDASSLOutput(const double &time,const double &h,const int &stp);

  ISolverSettings
    *_cppdasslsettings;              ///< Input      - Solver settings

  int
    _nrt,
    _idid,
    *_info,
    _lrw,
    _liw,
    *_iwork;
  double
    _atol,
    _rtol,
    _h,
    *_rwork,
    *_y,
    *_yp;
  int
    _dimSys;                 ///< Input       - (total) Dimension of system (=number of ODE)

  void
    *_data;

  static int res(const double* t, const double* y, const double* yprime, double* cj, double* delta, int* ires, void *par);
  int calcFunction(const double* t, const double* y, const double* yprime, double* cj, double* delta, int* ires);

  sparsematrix_t
    *_jac;

  // Variables for Coloured Jacobians
//  int  _sizeof_sparsePattern_colorCols;
//  int* _sparsePattern_colorCols;
//
//  int  _sizeof_sparsePattern_leadindex;
//  int* _sparsePattern_leadindex;
//
//
//  int  _sizeof_sparsePattern_index;
//  int* _sparsePattern_index;
//
//
//  int  _sparsePattern_maxColors;
//
//  bool _cvode_initialized;


//   ISystemProperties* _properties;
   IContinuous* _continuous_system[5];
//   IEvent* _event_system;
//   IMixedSystem* _mixed_system;
   ITime* _time_system[5];

//   std::vector<MeasureTimeData> measureTimeFunctionsArray;
//   MeasureTimeValues *measuredFunctionStartValues, *measuredFunctionEndValues;

};
#else
class CppDASSL : public ISolver, public SolverDefaultImplementation
{
public:
	CppDASSL(IMixedSystem* system, ISolverSettings* settings) : ISolver(), SolverDefaultImplementation(system, settings)
	{
		throw std::runtime_error("CppDASSL solver is not available.");
	}

	virtual void setStartTime(const double& time)
	{}

	virtual void setEndTime(const double& time)
	{}

	virtual void setInitStepSize(const double& stepSize)
	{}

	virtual void initialize()
	{}

	virtual bool stateSelection()
	{
		throw std::runtime_error("Peer solver is not available.");
	}

	virtual void solve(const SOLVERCALL command = UNDEF_CALL)
	{}

	virtual SOLVERSTATUS getSolverStatus()
	{ return UNDEF_STATUS; }

	virtual void setTimeOut(unsigned int time_out)
	{}

	virtual void stop()
	{}

	virtual void writeSimulationInfo()
	{}

};
#endif
