from multiprocessing import Process,Value,Manager
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
        
        self.result_list = Manager().list()
        """:obj:`list`: A list accessible from child processes; for the return values."""

        self.process_list= []
        """:obj:`list`: All of the child processes."""

        self.mem1 = Value("f",1)
        """:obj:`Value`: A value accessible from child processes; The memory needed to run a single child."""

        self.mem_a = Value("f",1)
        """:obj:`Value`: A value accessible from child processes; The current available memory."""

        self.mem_t = Value("f",psutil.virtual_memory()[1])
        """:obj:`Value`: A value accessible from child processes; The total amount of available memory at the creation of the instance."""

        self.nb_allowed = Value('i',1)
        """:obj:`Value`: The amount of theoretical child processes that can fit in FREE memory."""

        self.max_allowed = Value('i',1)
        """:obj:`Value`: The amount of theoretical child processes that can fit in TOTAL memory."""


        self.tasks = Value('i',0)
        """:obj:`Value`: A value accessible from child processes; The amount of child processes currently running."""

        self.started = Value('i',0)
        """:obj:`Value`: A value accessible from child processes; Number of started processes."""

        self.finished = Value('i',0)
        """:obj:`Value`: A value accessible from child processes; Number of finished processes."""

        self.amount = len(args)
        """:obj:`int`: The amount of processes that the instance will start."""
        
        self.start_time = 0
        """:obj:`float`: The present time."""

        self.time_for_one = Value("f",0)
        """:obj:`Value`: A value accessible from child processes; Time it takes for the first iteration to complete."""

        self.rewrite= {True:11,False:6}
        """:obj:`dict`: `int` value used in the terminal update rewriting process; amount of lines rewinded if large or small."""

    def execute_once(self,arg):
            """
            Method that is ran as a single .

            Parameters
            ----------
            args : :obj:`list`    A list of arguments given to the method

            Returns
            -------
            None    
            
            Note
            -----
            The return value of the passed method is passed to self.result_list
            """
            self.tasks.value+= 1
            self.started.value +=1
            self.processes.append(os.getpid())

            self.result_list.append(self.function(arg)) #execution of the passed method
            
            self.tasks.value -= 1
            self.finished.value += 1

    def processing_large_data(self):
        """
        Method for executing mutliprocess on tasks that demand a large amount of individual memory.
        Will, first, run a single process then, will start as many child processes as the available RAM permits
        Starting new ones as the RAM is freed

        Returns
        -------
        `list`
            The multiprocess-friendly list, that was updated by each child

        Notes
        -----
        If other application reduce the avalable RAM mid-execution,
        Multiprocess outputs "Killed" and kills the child.

        """
        print("    Starting multiprocessing, this might take some time\n    The first process is ran alone for calibration purposes")
        self.start_time = time.time()

        p = Process(target=self.execute_once, args=([self.args.pop(0)])) #Multiprocess runs once alone
        """:obj:`Process`: Representes a single child process"""

        p.start()

        print("\033[B"*self.rewrite[False], flush=True) #"\033[B" acts as a reverse \n, return the pointer one line up
        self.processes.append(os.getpid())              #adds the main thread on the list of processes to keep an eye on

        alfred = Process(target=self.buttler,args=([True]))  #ask the buttler to start complimentary processes
        alfred.start()
        
        time.sleep(1)   #give it a second to open the process, so it doesn't skip the while()
        while(self.tasks.value>0):
            #wait for the calibration process to finish
            time.sleep(0.1)
        self.time_for_one.value=time.time()-self.start_time
        
        while(len(self.args)!=0):
            if (self.tasks.value < self.max_allowed.value) & (self.nb_allowed.value >= 1):
                p = Process(target=self.execute_once, args=([self.args.pop(0)])) #Multiprocess runs the rest of the processes
                
                self.process_list.append(p)
                p.start()
                time.sleep(0.1)
            else:
                time.sleep(0.1)


        for p in self.process_list:
            p.join()
        time.sleep(1)   #give it a second to close the processes, do not remove

        alfred.terminate()
        print("\r", end="")
        print("\033[B"*self.rewrite[True], flush=True) #"\033[B" acts as a reverse \n, return the pointer one line up
        print("Completed with ",str(self.amount-self.finished.value)," errors\n") #if a process was killed or didn't finished; it will be know here

        return self.result_list

    def buttler(self,memBloc):
        """
        Ran as a child process, the buttler will constantly run other methods forever.
        In this case, it:
            - updates de memory capacity and 
            - prints updates on the terminal.
        It exists so not to bottleneck the main thread.

        Note
        ----
        Uses timers to execute it's methods because time.sleep() it processor hungry if constantly called
        """
        terminal = time.time()
        """:obj:`float`: Time display for 0.1s lapse."""
        mem = time.time()
        """:obj:`float`: Time display for 1s lapse."""
        while True:
            now = time.time()
            if now-terminal>0.1:    #has 0.1 second passed?
                self.terminal_update(memBloc)
                terminal = now
            if now-mem>1 & memBloc: #has 1 second passed?
                self.mem_update()
                mem = now

    def mem_update(self):
        """
        Method that sets the baseline for memory calculation
        + output some information to the terminal
        All memory values are in bytes

        This method is ran from the `buttler()` and updates every second

        """

        memBuffer = 0.9 #%
        """:obj:`float`: %Amount of bytes to substract from the available RAM for safety purposes."""
        self.mem_a.value = psutil.virtual_memory()[1]*memBuffer
        for child in self.processes:
            try:    #in a try/except because processes ID are never removed from the list
                mem = psutil.Process(child).memory_full_info()[8] #uss memory usage; humch much is this process using NOW
                """:obj:`float`: Amount of bytes."""
                if (self.mem1.value<mem):
                    self.mem1.value=mem #does it for the whole run in case this maximum is increased by future childs
            except:
                ""

        self.nb_allowed.value = math.floor(((self.mem_a.value/self.mem1.value)))
        if self.nb_allowed.value < 1:
            self.nb_allowed.value = 1    #Need to at least be able to start a single process
        self.max_allowed.value = math.floor(((self.mem_t.value/self.mem1.value)))
        if self.max_allowed.value < 1:
            self.max_allowed.value = 1   #Need to at least be able to start a single process

    def terminal_update(self, memBlock):
        """
        Method that constantly updates the user about the currently run tasks

        Note
        ----
        This method is ran from the buttler() and updates every 0.1 second
        """

        s=self.started.value
        a=self.amount
        f=self.finished.value
        nowTime = time.time()
        eTime = round((nowTime - self.start_time)*10)/10 #current execution time
        
        print("---")
        if(memBlock): #block of prints used only by processing_large_data()
            print("Available memory: ", round(self.mem_a.value/10000000)/100,"/",round(self.mem_t.value/10000000)/100, "Gb        ", end="\n", flush=True)
            print("Active processes: ", str(self.tasks.value) + " / " + str(self.max_allowed.value) + "            ", end="\n", flush=True)
            print("Min memory per:   ",round(self.mem1.value/10000000)/100, "Gb        ", end="\n", flush=True)
            print("Time for one:     " ,round(self.time_for_one.value*10)/10," seconds               ", flush=True)
            print("---")
        print("Started:          ",s,"/",a,"   ",round(s/a*100),"%           ", flush=True)
        print("Finished:         ",f,"/",s,"   ",round(f/a*100),"%           ", flush=True)
        print("Time elapsed:     " ,eTime," seconds               ", flush=True)
        print("---")

        print("\r", end="",flush=True)
        print("\033[A"*self.rewrite[memBlock], flush=True) #"\033[B" acts as a reverse \n, return the pointer one line up

    def execute_small(self,arg):
        """
        The method executed by processingSmallData
        
        Note:
        --------
        No returns, but the return value of the executed method is passed to a global multiprocessing-friendly list
        """
        self.started.value +=1

        result = self.function(arg)
        self.result_list.append(result)

        self.finished.value += 1
    
    def processing_small_data(self):
        """
        Method for executing mutliprocess on tasks that demand little to no memory.
        Will immedialty start all the child processes.
        The packaging causes some marginal time lost; 
        Only use for methods that take at least a second to run
        below that, a for loop is likely much faster.
        

        Variables:
            p   Process   Representes a single child process
            a   None      Exists only to permit the for loop

        Returns:
        result_list : :obj:`list`
            The multiprocess-friendly list, that was updated by each child

        Notes
        ----
        Errors if other application reduce the avalable RAM mid-execution.
        Multiprocess outputs "Killed" and kills the child.

        """

        self.start_time = time.time()
        alfred = Process(target=self.buttler,args=([False]))
        alfred.start()

    
        for a in range(len(self.args)):
            p = Process(target=self.execute_small, args=([self.args.pop(0)])) #Multiprocess runs
            self.process_list.append(p)
            p.start()

        for p in self.process_list:
            p.join()

        time.sleep(1) #wait for the processes to close; do not remove
        finishedTime = round(time.time()-self.start_time)*10/10

        alfred.terminate()
        print("\r", end="")
        print("\033[B"*self.rewrite[False], flush=True) #"\033[B" acts as a reverse \n, return the pointer one line up
        print("Completed ",len(self.result_list),"tasks in ",finishedTime,"seconds")
 
        return self.result_list
