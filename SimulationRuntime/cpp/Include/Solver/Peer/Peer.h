
#pragma once

#include "FactoryExport.h"

#include <Core/Solver/SolverDefaultImplementation.h>
#include <Core/Utils/extension/measure_time.hpp>


/*****************************************************************************/
// Peer
// BDF-Verfahren für steife und nicht-steife ODEs
// Dokumentation siehe offizielle Peer Doku

/*****************************************************************************
Copyright (c) 2014, IWR TU Dresden, All rights reserved
*****************************************************************************/
class Peer
  : public ISolver,  public SolverDefaultImplementation
{
public:

  Peer(IMixedSystem* system, ISolverSettings* settings);

  virtual ~Peer();

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
  void writePeerOutput(const double &time,const double &h,const int &stp);
  void evalJ(const double& t, const double* z, double* T, double fac=1);
  void evalF(const double& t, const double* z, double* f);
  void evalD(const double& t, const double* y, double* T);
  void setcycletime(double cycletime);
  void ros2(double * y, double& tstart, double tend);

  ISolverSettings
    *_peersettings;              ///< Input      - Solver settings

  bool
    _first;

  long int
    _dimSys;                 ///< Input       - (total) Dimension of system (=number of ODE)


    int
        _rstages,
        _rank,
        _size;

    long int
        *_P;

    double
        *_G,
        *_E,
        *_Theta,
        *_c,
        *_F,
        *_y,
        *_Y1,
        *_Y2,
        *_Y3,
        *_T,
        _h;



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
   IContinuous* _continuous_system;
//   IEvent* _event_system;
//   IMixedSystem* _mixed_system;
   ITime* _time_system;

//   std::vector<MeasureTimeData> measureTimeFunctionsArray;
//   MeasureTimeValues *measuredFunctionStartValues, *measuredFunctionEndValues;

};

