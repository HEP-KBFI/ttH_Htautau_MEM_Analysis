# include "ttH_Htautau_MEM_Analysis/MEM/interface/PyConfig.h"
# include <algorithm>
# include <iostream> 

using namespace std;

static void printPyConfigUsage() {
   cerr << "./xxxx.exe .... py_config_file_name[.py]" << endl;
}

string fileToPyModuleName( const char *filename ) {
  string name( filename ); 
  size_t pos = name.find(".py", name.length()-3);
  if ( pos != string::npos ) {
    name.resize( name.length()-3 );    
  }
  // '/' -> '.'
  std::replace( name.begin(), name.end(), '/', '.');
  // cerr << "file name: " << name << endl;
  
  return name;
}

void debugPyConfigFile(PyObject* pModule ) {
  PyObject *pDict = 0;
  int i;
  
  if (pModule != NULL) {
    pDict = PyModule_GetDict(pModule);
    PyObject *pKeys = PyDict_Keys(pDict); 
    // Get the variable or identifier in python file
    for( i = 0; i < PyList_Size(pKeys); ++i) {
      PyObject *pKey = PyList_GetItem(pKeys, i); 
      PyObject *pValue = PyDict_GetItem(pDict, pKey); 
      // Debug 
      cout << "KEY "<< PyString_AsString(pKey) << endl;
      assert(pValue);
    }
    Py_DECREF(pKeys);
  }
  Py_DECREF(pDict);

}

PyObject* readPyConfig( const char* filename ) {

  PyObject *pName, *pModule;
  string fileName = fileToPyModuleName( filename );

  Py_Initialize();
  if (fileName.size() != 0 ) 
    pName = PyString_FromString( fileName.c_str() );
  else 
    pName = 0;

  /* Error checking of pName */
  if( pName == 0 ) {
    cerr << "PyConfig: Bad configuration file name:" 
         << fileName << endl;
    printPyConfigUsage();
    return 0;
  } 

  pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  
  if( pModule == 0 ) {
    cerr << "PyConfig: Bad python module name " 
         << fileName << endl;
    printPyConfigUsage();
    return 0;
  } 
  return pModule;
}

 
void closePyConfig( PyObject* pModule ) {

  Py_DECREF(pModule);
  Py_Finalize();
}

double getPyObjetDouble( PyObject* pyObj, const char *name) {
  double value=nan("NaN");
  PyObject* pVariable = PyObject_GetAttrString( pyObj, name);
  if (pVariable != 0 ) {
    value = PyFloat_AsDouble( pVariable );
  }
  return value;
}

double getPyDictDouble( PyObject* pyDict, const char *name) {
  double value=nan("NaN");
  PyObject* pItem = PyDict_GetItemString(pyDict, name);
  if (pItem != 0 ) {
    value = PyFloat_AsDouble( pItem );
  }
  return value;
}

float getPyDictFloat( PyObject* pyDict, const char *name) {
  float value=nan("NaN");
  PyObject* pItem = PyDict_GetItemString(pyDict, name);
  if (pItem != 0 ) {
    value = (float) (PyFloat_AsDouble( pItem ));
  }
  return value;
}

long int getPyDictLong( PyObject* pyDict, const char *name) {
  long int value = LONG_MAX;
  PyObject* pItem = PyDict_GetItemString(pyDict, name);
  if (pItem != 0 ) {
    value = PyLong_AsLong( pItem );
  }
  return value;
}

int getPyDictInt( PyObject* pyDict, const char *name) {
  int value = INT_MAX;
  PyObject* pItem = PyDict_GetItemString(pyDict, name);
  if (pItem != 0 ) {
    value = PyLong_AsLong( pItem );
  }
  return value;
}



vector<float> getPyDictVectFloat(PyObject* pyDict, const char *name) {

  vector<float> vect;

  PyObject* pVariable = PyDict_GetItemString(pyDict, name);
  if (pVariable!=0){
    for(int i=0; i<PyList_Size(pVariable); i++){
      PyObject* pItem = PyList_GetItem(pVariable,i);    
      vect.push_back((float) PyFloat_AsDouble(pItem));
    }
  }
  return vect;

}




vector<int> getPyDictVectInt(PyObject* pyDict, const char *name) {

  vector<int> vect;

  PyObject* pVariable = PyDict_GetItemString(pyDict, name);
  if (pVariable!=0){
    for(int i=0; i<PyList_Size(pVariable); i++){
      PyObject* pItem = PyList_GetItem(pVariable,i);
      vect.push_back((int) PyLong_AsLong(pItem));      
    }
  }

  return vect;

}
