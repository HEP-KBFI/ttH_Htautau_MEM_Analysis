/* 
 * File:   ConfigFile.h
 * Author: grasseau
 *
 * Created on 25 novembre 2014, 17:57
 */

#ifndef PYCONFIG_H
#define	PYCONFIG_H

#include <boost/python.hpp>
#include <vector>


PyObject* readPyConfig     ( const char* filename );
void      debugPyConfigFile( PyObject* pModule );
void      closePyConfig    ( PyObject* pModule );

double getPyObjetDouble( PyObject* pyObj, const char *name);

double getPyDictDouble( PyObject* pyDict, const char *name);

float getPyDictFloat( PyObject* pyDict, const char *name);

long int getPyDictLong( PyObject* pyDict, const char *name);

int getPyDictInt( PyObject* pyDict, const char *name);

std::vector<float> getPyDictVectFloat( PyObject* pyDict, const char *name);

std::vector<int> getPyDictVectInt( PyObject* pyDict, const char *name);
#endif	/* PYCONFIG_H */

