/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-2014, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF THE BSD NEW LICENSE OR THE
 * GPL VERSION 3 LICENSE OR THE OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the OSMC (Open Source Modelica Consortium)
 * Public License (OSMC-PL) are obtained from OSMC, either from the above
 * address, from the URLs: http://www.openmodelica.org or
 * http://www.ida.liu.se/projects/OpenModelica, and in the OpenModelica
 * distribution. GNU version 3 is obtained from:
 * http://www.gnu.org/copyleft/gpl.html. The New BSD License is obtained from:
 * http://www.opensource.org/licenses/BSD-3-Clause.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, EXCEPT AS
 * EXPRESSLY SET FORTH IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE
 * CONDITIONS OF OSMC-PL.
 *
 */

#include "stateset.h"
#include "util/omc_error.h"

#include <memory.h>

/*! \fn initializeStateSetJacobians
 *
 *  initialize jacobians for state selection
 *
 *  \param [ref] [data] ???
 *
 *  \author ???
 */
void initializeStateSetJacobians(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  long i = 0;
  STATE_SET_DATA *set = NULL;

  /* go troug all state sets*/
  for(i=0; i<data->modelData->nStateSets; i++)
  {
    set = &(data->simulationInfo->stateSetData[i]);
    if(set->initialAnalyticalJacobian(data, threadData))
    {
      throwStreamPrint(threadData, "can not initialze Jacobians for dynamic state selection");
    }
  }
  initializeStateSetPivoting(data);
  TRACE_POP
}

/*! \fn initializeStateSetPivoting
 *
 *  initialize pivoting data for state selection
 *
 *  \param [ref] [data] ???
 *
 *  \author ???
 */
void initializeStateSetPivoting(DATA *data)
{
  TRACE_PUSH
  long i = 0;
  long n = 0;
  STATE_SET_DATA *set = NULL;
  unsigned int aid = 0;
  modelica_integer *A = NULL;

  /* go troug all state sets */
  for(i=0; i<data->modelData->nStateSets; i++)
  {
    set = &(data->simulationInfo->stateSetData[i]);
    aid = set->A->id - data->modelData->integerVarsData[0].info.id;
    A = &(data->localData[0]->integerVars[aid]);

    memset(A, 0, set->nCandidates*set->nStates*sizeof(modelica_integer));

    /* initialize row and col indizes */
    for(n=0; n<set->nDummyStates; n++)
      set->rowPivot[n] = n;

    for(n=0; n<set->nCandidates; n++)
      set->colPivot[n] = set->nCandidates-n-1;

    for(n=0; n<set->nStates; n++)
      A[n*set->nCandidates + n] = 1;  /* set A[row, col] */
  }
  TRACE_POP
}

/*! \fn freeStateSetData
 *
 *  free jacobians for state selection
 *
 *  \param [ref] [data] ???
 *
 *  \author ???
 */
void freeStateSetData(DATA *data)
{
  TRACE_PUSH
  long i=0;

  /* go through all state sets */
  for(i=0; i<data->modelData->nStateSets; i++)
  {
     STATE_SET_DATA *set = &(data->simulationInfo->stateSetData[i]);
     free(set->states);
     free(set->statescandidates);
     free(set->rowPivot);
     free(set->colPivot);
     free(set->J);
  }
  TRACE_POP
}

/*! \fn getAnalyticalJacobianSet
 *
 *  function calculates analytical jacobian
 *
 *  \param [ref] [data] ???
 *  \param [out] [index] ???
 *
 *  \author wbraun
 */
static void getAnalyticalJacobianSet(DATA* data, threadData_t *threadData, unsigned int index)
{
  TRACE_PUSH
  unsigned int i, j, k, l, ii;
  unsigned int jacIndex = data->simulationInfo->stateSetData[index].jacobianIndex;
  unsigned int nrows = data->simulationInfo->analyticJacobians[jacIndex].sizeRows;
  unsigned int ncols = data->simulationInfo->analyticJacobians[jacIndex].sizeCols;
  double* jac = data->simulationInfo->stateSetData[index].J;

  /* set all elements to zero */
  memset(jac, 0, (nrows*ncols*sizeof(double)));

  for(i=0; i < data->simulationInfo->analyticJacobians[jacIndex].sparsePattern.maxColors; i++)
  {
    for(ii=0; ii < data->simulationInfo->analyticJacobians[jacIndex].sizeCols; ii++)
      if(data->simulationInfo->analyticJacobians[jacIndex].sparsePattern.colorCols[ii]-1 == i)
        data->simulationInfo->analyticJacobians[jacIndex].seedVars[ii] = 1;

    if(ACTIVE_STREAM(LOG_DSS_JAC))
    {
      infoStreamPrint(LOG_DSS_JAC, 1, "Caluculate one col:");
      for(l=0; l < data->simulationInfo->analyticJacobians[jacIndex].sizeCols; l++)
        infoStreamPrint(LOG_DSS_JAC, 0, "seed: data->simulationInfo->analyticJacobians[index].seedVars[%d]= %f", l, data->simulationInfo->analyticJacobians[jacIndex].seedVars[l]);
      messageClose(LOG_DSS_JAC);
    }

    (data->simulationInfo->stateSetData[index].analyticalJacobianColumn)(data, threadData);

    for(j=0; j < data->simulationInfo->analyticJacobians[jacIndex].sizeCols; j++)
    {
      if(data->simulationInfo->analyticJacobians[jacIndex].seedVars[j] == 1)
      {
        if(j==0)
          ii = 0;
        else
          ii = data->simulationInfo->analyticJacobians[jacIndex].sparsePattern.leadindex[j-1];

        /* infoStreamPrint(LOG_DSS_JAC, 0, "take for %d -> %d\n", j, ii); */

        while(ii < data->simulationInfo->analyticJacobians[jacIndex].sparsePattern.leadindex[j])
        {
          l  = data->simulationInfo->analyticJacobians[jacIndex].sparsePattern.index[ii];
          k  = j*data->simulationInfo->analyticJacobians[jacIndex].sizeRows + l;
          jac[k] = data->simulationInfo->analyticJacobians[jacIndex].resultVars[l];
          /* infoStreamPrint(LOG_DSS_JAC, 0, "write %d. in jac[%d]-[%d, %d]=%f from col[%d]=%f", ii, k, l, j, jac[k], l, data->simulationInfo->analyticJacobians[jacIndex].resultVars[l]); */
          ii++;
        };
      }
    }
    for(ii=0; ii < data->simulationInfo->analyticJacobians[jacIndex].sizeCols; ii++)
      if(data->simulationInfo->analyticJacobians[jacIndex].sparsePattern.colorCols[ii]-1 == i)
        data->simulationInfo->analyticJacobians[jacIndex].seedVars[ii] = 0;
  }

  /*
  if(ACTIVE_STREAM(LOG_DSS))
  {
    char buffer[4096];

    infoStreamPrint(LOG_DSS, "jacobian %dx%d [id: %d]", data->simulationInfo->analyticJacobians[jacIndex].sizeRows, data->simulationInfo->analyticJacobians[jacIndex].sizeCols, jacIndex);
    INDENT(LOG_DSS);
    for(i=0; i<data->simulationInfo->analyticJacobians[jacIndex].sizeRows; i++)
    {
      buffer[0] = 0;
      for(j=0; j < data->simulationInfo->analyticJacobians[jacIndex].sizeCols; j++)
        sprintf(buffer, "%s%.5e ", buffer, jac[i*data->simulationInfo->analyticJacobians[jacIndex].sizeCols+j]);
      infoStreamPrint(LOG_DSS, "%s", buffer);
    }
    RELEASE(LOG_DSS);
  }
  */
  TRACE_POP
}

/*! \fn setAMatrix
 *
 *  ??? desc ???
 *
 *  \param [ref] [newEnable]
 *  \param [ref] [nCandidates]
 *  \param [ref] [nStates]
 *  \param [ref] [Ainfo]
 *  \param [ref] [states]
 *  \param [ref] [statecandidates]
 *  \param [ref] [data]
 */
static void setAMatrix(modelica_integer* newEnable, modelica_integer nCandidates, modelica_integer nStates, VAR_INFO* Ainfo, VAR_INFO** states, VAR_INFO** statecandidates, DATA *data)
{
  TRACE_PUSH
  modelica_integer col;
  modelica_integer row=0;
  /* clear old values */
  unsigned int aid = Ainfo->id - data->modelData->integerVarsData[0].info.id;
  modelica_integer *A = &(data->localData[0]->integerVars[aid]);
  memset(A, 0, nCandidates*nStates*sizeof(modelica_integer));

  for(col=0; col<nCandidates; col++)
  {
    if(newEnable[col]==2)
    {
      unsigned int firstrealid = data->modelData->realVarsData[0].info.id;
      unsigned int id = statecandidates[col]->id-firstrealid;
      unsigned int sid = states[row]->id-firstrealid;
      infoStreamPrint(LOG_DSS, 0, "select %s", statecandidates[col]->name);
      /* set A[row, col] */
      A[row*nCandidates + col] = 1;
      /* reinit state */
      data->localData[0]->realVars[sid] = data->localData[0]->realVars[id];
      row++;
    }
  }
  TRACE_POP
}

/*! \fn comparePivot
 *
 *  ??? desc ???
 *
 *  \param [ref] [oldPivot]
 *  \param [ref] [newPivot]
 *  \param [ref] [nCandidates]
 *  \param [ref] [nDummyStates]
 *  \param [ref] [nStates]
 *  \param [ref] [A]
 *  \param [ref] [states]
 *  \param [ref] [statecandidates]
 *  \param [ref] [data]
 *  \return ???
 */
static int comparePivot(modelica_integer *oldPivot, modelica_integer *newPivot, modelica_integer nCandidates, modelica_integer nDummyStates, modelica_integer nStates, VAR_INFO* A, VAR_INFO** states, VAR_INFO** statecandidates, DATA *data, int switchStates)
{
  TRACE_PUSH
  modelica_integer i;
  int ret = 0;
  modelica_integer* oldEnable = (modelica_integer*) calloc(nCandidates, sizeof(modelica_integer));
  modelica_integer* newEnable = (modelica_integer*) calloc(nCandidates, sizeof(modelica_integer));

  for(i=0; i<nCandidates; i++)
  {
    modelica_integer entry = (i < nDummyStates) ? 1: 2;
    newEnable[ newPivot[i] ] = entry;
    oldEnable[ oldPivot[i] ] = entry;
 }

  for(i=0; i<nCandidates; i++)
  {
    if(newEnable[i] != oldEnable[i])
    {
      if(switchStates)
      {
        infoStreamPrint(LOG_DSS, 1, "select new states at time %f", data->localData[0]->timeValue);
        setAMatrix(newEnable, nCandidates, nStates, A, states, statecandidates, data);
        messageClose(LOG_DSS);
      }
      ret = -1;
      break;
    }
  }

  free(oldEnable);
  free(newEnable);

  TRACE_POP
  return ret;
}

/*! \fn stateSelection
 *
 *  function to select the actual states
 *
 *  \param [ref] [data]
 *  \param [in]  [reportError]
 *  \param [in]  [switchStates] flag for switch states, function does switch only if this switchStates = 1
 *  \return ???
 *
 *  \author Frenkel TUD
 */
int stateSelection(DATA *data, threadData_t *threadData, char reportError, int switchStates)
{
  TRACE_PUSH
  long i=0;
  long j=0;
  int globalres=0;

  /* go through all the state sets */
  for(i=0; i<data->modelData->nStateSets; i++)
  {
    int res=0;
    STATE_SET_DATA *set = &(data->simulationInfo->stateSetData[i]);
    modelica_integer* oldColPivot = (modelica_integer*) malloc(set->nCandidates * sizeof(modelica_integer));
    modelica_integer* oldRowPivot = (modelica_integer*) malloc(set->nDummyStates * sizeof(modelica_integer));

    /* generate jacobian, stored in set->J */
    getAnalyticalJacobianSet(data, threadData, i);

    /* call pivoting function to select the states */
    memcpy(oldColPivot, set->colPivot, set->nCandidates*sizeof(modelica_integer));
    memcpy(oldRowPivot, set->rowPivot, set->nDummyStates*sizeof(modelica_integer));
    if((pivot(set->J, set->nDummyStates, set->nCandidates, set->rowPivot, set->colPivot) != 0) && reportError)
    {
      /* error, report the matrix and the time */
      char buffer[4096];

      warningStreamPrint(LOG_DSS, 1, "jacobian %dx%d [id: %ld]", data->simulationInfo->analyticJacobians[set->jacobianIndex].sizeRows, data->simulationInfo->analyticJacobians[set->jacobianIndex].sizeCols, set->jacobianIndex);
      for(i=0; i < data->simulationInfo->analyticJacobians[set->jacobianIndex].sizeRows; i++)
      {
        buffer[0] = 0;
        for(j=0; j < data->simulationInfo->analyticJacobians[set->jacobianIndex].sizeCols; j++)
          sprintf(buffer, "%s%.5e ", buffer, set->J[i*data->simulationInfo->analyticJacobians[set->jacobianIndex].sizeCols+j]);
        warningStreamPrint(LOG_DSS, 0, "%s", buffer);
      }

      for(i=0; i<set->nCandidates; i++)
        warningStreamPrint(LOG_DSS, 0, "%s", set->statescandidates[i]->name);
      messageClose(LOG_DSS);

      throwStreamPrint(threadData, "Error, singular Jacobian for dynamic state selection at time %f\nUse -lv LOG_DSS_JAC to get the Jacobian", data->localData[0]->timeValue);
    }
    /* if we have a new set throw event for reinitialization
       and set the A matrix for set.x=A*(states) */
    res = comparePivot(oldColPivot, set->colPivot, set->nCandidates, set->nDummyStates, set->nStates, set->A, set->states, set->statescandidates, data, switchStates);
    if(!switchStates)
    {
      memcpy(set->colPivot, oldColPivot, set->nCandidates*sizeof(modelica_integer));
      memcpy(set->rowPivot, oldRowPivot, set->nDummyStates*sizeof(modelica_integer));
    }
    if(res)
      globalres = 1;

    free(oldColPivot);
    free(oldRowPivot);
  }

  TRACE_POP
  return globalres;
}
