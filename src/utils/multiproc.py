from multiprocess import Process,Value,Manager
import psutil
import time
import math
import os


class Multi:
    """

    Class that harnesses the power of "multiprocess" and wraps it as a single callable line.
    It allows to easily run RAM-heavy and multiple tasks on servers.

    """

    def __init__(self,args,function):
        """
        Constructor of the multiprocessing tool 

        Attributes
        ----------
        args : :obj:`list`
            The list of value to run. One process will be born per 1st level member
        function : :obj:`function`
            The `function` that needs to be run as a multiple instance
        
        Note
        ----
        The constructor needs two thing, a `function`, and a `list` of `lists`. 
        Each `sub-list` of the `list` will be the `list` of arguments given to the child process running the given `function`
        **The function used must only take a single list as input; and then, deal with it's division into other variables**
        
        Example
        -------
        ex.: Run 3 processes of g()

            def g(args):
                n1=args[0]
                n2=args[1]
                return n1+n2
            
            list = [ [0,1], [1,1] ,[2,1] ]

            mps = Multi( list, g ).processingSmallData()        #note that "g" is not written "g()" we want to pass it, not run it
            print( mps )

            >>[1,2,3]

        """

        self.args = Manager().list(args)
        self.function = function
        
        self.processes = Manager().list()
        """:obj:`list`: A list accessible from child processes; list of all the process ID started in this object."""
        
        self.resultList = Manager().list()
        """:obj:`list`: A list accessible from child processes; for the return values."""

        self.processlist= []
        """:obj:`list`: All of the child processes."""

        self.mem1 = Value("f",1)
        """:obj:`ctypes`: A value accessible from child processes; The memory needed to run a single child."""

        self.memA = Value("f",1)
        """:obj:`ctypes`: A value accessible from child processes; The current available memory."""

        self.memT = Value("f",psutil.virtual_memory()[1])
        """:obj:`ctypes`: A value accessible from child processes; The total amount of available memory at the creation of the instance."""

        self.nbAllowed = Value('i',1)
        """:obj:`ctypes`: The amount of theoretical child processes that can fit in FREE memory."""

        self.maxAllowed = Value('i',1)
        """:obj:`ctypes`: The amount of theoretical child processes that can fit in TOTAL memory."""


        self.tasks = Value('i',0)
        """:obj:`ctypes`: A value accessible from child processes; The amount of child processes currently running."""

        self.started = Value('i',0)
        """:obj:`ctypes`: A value accessible from child processes; Number of started processes."""

        self.finished = Value('i',0)
        """:obj:`ctypes`: A value accessible from child processes; Number of finished processes."""

        self.amount = len(args)
        """:obj:`int`: The amount of processes that the instance will start."""
        
        self.startTime = 0
        """:obj:`float`: The present time."""

        self.timeForOne = Value("f",0)
        """:obj:`ctypes`: A value accessible from child processes; Time it takes for the first iteration to complete."""

        self.rewrite= {True:11,False:6}
        """:obj:`dict`: `int` value used in the terminal update rewriting process; amount of lines rewinded if large or small."""
