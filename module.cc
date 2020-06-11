#include <iostream>
#include <sys/resource.h> 
#include <sys/time.h> 
#include <unistd.h> 
#include <stdio.h> 
#include <time.h>
#include<bits/stdc++.h>
#include <sys/stat.h> 
#include <sys/types.h> 
#define PY_SSIZE_T_CLEAN
#include <Python/Python.h>
#include <pyhelper.hpp>

vector<int> listTupleToVector_Int(CPyObject* incoming) {
	vector<int> data;
	if (PyTuple_Check(incoming)) {
		for(Py_ssize_t i = 0; i < PyTuple_Size(incoming); i++) {
			PyObject *value = PyTuple_GetItem(incoming, i);
			data.push_back( PyFloat_AsDouble(value) );
		}
	} else {
		if (PyList_Check(incoming)) {
			for(Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
				PyObject *value = PyList_GetItem(incoming, i);
				data.push_back( PyFloat_AsDouble(value) );
			}
		} else {
			throw logic_error("Passed PyObject pointer was not a list or tuple!");
		}
	}
	return data;
}

int main(int argc, char* argv[]) {
    CPyInstance hInstance;

	CPyObject pName = PyUnicode_FromString("pyemb3");
	CPyObject pModule = PyImport_Import(pName);

	if(pModule)
	{
		CPyObject pFunc = PyObject_GetAttrString(pModule, "getInteger");
		if(pFunc && PyCallable_Check(pFunc))
		{
			CPyObject pValue = PyObject_CallObject(pFunc, NULL);

			vector<int> resVec = listTupleToVector_Int(pValue)
            for (auto it : resvec) {
                std::cout << *it << std::endl;
            }
		}
		else
		{
			printf("ERROR: function getInteger()\n");
		}

	}
	else
	{
		printf_s("ERROR: Module not imported\n");
	}

	return 0;

}