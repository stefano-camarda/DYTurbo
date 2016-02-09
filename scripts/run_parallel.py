#!/opt/rh/python27/root/usr/bin/python
# -*- coding: utf-8 -*- 

## Documentation for file
#
# More details. 
#
# @file run_paralel.py
# @author cuto <Jakub.Cuth@cern.ch>
# @date 2015-04-16


import os, shutil, sys, errno, time, threading


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
        pass
    pass

DRYRUN=True
DRYRUN=False
class BashCommand(threading.Thread) :
    def __init__(s, _command="echo yes", _jname="bash", _tmstmp="",_seed=66) :
        s.command=_command
        timestamp=time.strftime("%s") if _tmstmp=="" else _tmstmp
        s.seed=_seed
        #basefolder = os.path.expanduser("~")
        basefolder = os.getcwd()
        s.logDir_base=basefolder+"/parallel/log"
        s.runDir_base=basefolder+"/parallel/run"
        s.jobName="{}_{}_{}_{}".format(_jname,os.environ["USER"],s.seed,timestamp )
        # define names
        threading.Thread.__init__(s,name=s.jobName)
        s.runDir=s.runDir_base+"/"+s.jobName
        s.logStd=s.logDir_base+"/"+s.jobName+".out"
        s.logErr=s.logDir_base+"/"+s.jobName+".err"
        s.rc=-999
        pass

    def pr(s,text):
        print "* {} : {}".format(s.jobName,text)
        pass

    def init(s):
        # check dirs
        mkdir_p(s.logDir_base)
        mkdir_p(s.runDir)
        # init log
        f_out = open(s.logStd,'w')
        f_out.write("JobName:     {}\n"         .format( s.jobName                              ))
        f_out.write("System info: {}\n"         .format(  " ".join(os.uname())                  ))
        f_out.write("Start time:  {}\n"         .format( time.strftime("%a, %d %b %Y %H:%M:%S") ))
        f_out.write("Command: \n```\n{}\n```\n" .format( s.command                              ))
        f_out.write("\n\n" + 5*"-" + "JOBSTART" + 5*"-" + "\n")
        pass

    def finalize(s):
        # final log
        f_out = open(s.logStd,'a')
        f_out.write("\n\n" + 5*"-" + "JOBEND" + 5*"-" + "\n")
        f_out.write("Return code: {}\n" .format( s.jobName                              ))
        f_out.write("End time:    {}\n" .format( time.strftime("%a, %d %b %Y %H:%M:%S") ))
        # rm dirs
        #os.removedirs(s.runDir)
        shutil.rmtree(s.runDir)
        pass

    def run(s):
        s.init()
        # run
        s.pr("job start")
        query="cd {}; export LSB_JOBINDEX={}; {} >> {} 2>> {}".format( s.runDir, s.seed, s.command, s.logStd, s.logErr)
        if DRYRUN :
            print query
        else :
            s.rc = os.system(query)
        s.finalize()
        s.pr("job end")
        pass

def run_all_command( command_list, jobname_list, seedslist, N_threads=15):
    subm_time= time.strftime("%s")
    for i_job,comm in enumerate(command_list) :
        jname = jobname_list[i_job]
        seedsStr = seedslist[i_job]
        seeds = [ 100  ]
        if len(seedsStr)==0:
            pass
        elif "-" in seedsStr :
            seeds=range(
                    int(seedsStr.split("-")[0]),
                    int(seedsStr.split("-")[1]),
                    )
        elif "," in seedsStr :
            seeds= [int(x) for x in seedsStr.split(",")]
        else :
            seeds= [int(seedsStr)]
        # wait until some slot is empty (active jobs < Nthreads )
        for seed in seeds :
            time.sleep(1)
            while threading.active_count() > N_threads :
                time.sleep(1)
            fjob = BashCommand(comm,jname,subm_time,seed)
            fjob.start()
            pass
        pass
    pass

def readinput(fname="") :
    commands=list()
    names=list()
    seedlist=list()
    cmdlist=sys.stdin
    if fname!="" : cmdlist=open(fname)
    for line in cmdlist:
        stripped = line.strip()
        if not stripped: break
        parse_space = stripped.split(" ")
        if len(parse_space) == 1 :
            names    .append("bash")
            commands .append(parse_space[0])
        elif len(parse_space) == 2 :
            names    .append(parse_space[0])
            commands .append(parse_space[1])
            seedlist .append("123456")
        elif len(parse_space) == 3 :
            names    .append(parse_space[0])
            commands .append(parse_space[1])
            seedlist .append(parse_space[2])
            pass
        pass
    #print commands
    #print names
    #print seedlist
    # todo seeds
    return commands,names,seedlist

## Documentation for main
#
# More details. 
if __name__ == '__main__' :
    commands,names,seedslist = readinput("scripts/cmd_list")
    ##testing = [ "echo test {} ".format(x) for x in range(0,30)]
    run_all_command(commands, names, seedslist, 3)
    pass


