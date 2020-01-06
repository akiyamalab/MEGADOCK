/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Mpidp
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Mpidp_h
#define Mpidp_h 1

#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <mpi.h>

using namespace std;

// JOB management table
typedef struct {
    string	    name;		// Job name
    int		    exec;		// EXEC
    int		    status;		// calculation control flag
    vector<int>	worker;		// WID
    vector<int>	rcode[3];	// 0:END 1:RET 2:FILE
} NameLog;

// Worker management table
typedef struct {
    vector<string> name;		// Job name
    vector<int> 	rcode;		// 0: failure 1:success
    int		failure;	// calculation failure counter
} WorkerLog;

class Mpidp
{
 private:
    Mpidp(Mpidp &c){}
    const Mpidp & operator=(const Mpidp &c);
    MPI_Status      _Status;
    string          _Table_file;
    string          _Out_file;
    int             _Out_option;
    int             _Ntry;
    int			    _Worker_life;
    string		    _Title;
    string		    _Param;
    char			_Name[7];
    vector<string>	_Table_list;
    NameLog		    *_Namelog;
    WorkerLog       *_Workerlog;
    int			    _Csize;
    int			    _Psize;
    int			    _Ndata;
 protected:
    virtual string	erase_space(const string &s0,const int ip);
    virtual int		argument(int argc,char *argv[],char **wargv);
    virtual int		argument(int argc,char *argv[],string &main_argv);
    virtual int		for_worker(int &retry,char *ctable,int argc2,
                               char **wargv);
    virtual void	for_worker(int &retry,char *ctable,int &ia,
                               string &argv_joblist);
    virtual string	replace_pattern(const string &pattern,
                                    const string &position,
                                    const string &option);
 public:
    Mpidp() {
#ifdef DEBUG
        cout << "Constructing Mpidp.\n";
#endif
    }
    virtual ~Mpidp() {
#ifdef DEBUG
        cout << "Destructing Mpidp.\n";
#endif
    }
    virtual void read_table(int argc,char *argv[],int &ntry,
                            ofstream &logout);
    virtual int  master0(const int &nproc);
    virtual int	 master(const int &nproc);
    virtual void worker(int &myid,char *hostname,int argc,
                        char *argv[]);
    virtual void write_table(const int &nproc,ofstream &logout);
};

#endif
