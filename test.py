#!/usr/bin/env python3

import multiprocessing
import time
import random
import os
from multiprocessing import Queue

# q = Queue()

def hello(n, q):
    time.sleep(random.randint(1,3))
    q.put(os.getpid())
    print("[{0}] Hello!".format(n))

def go():
    q = Queue()
    processes = [ ]
    for i in range(10):
        t = multiprocessing.Process(target=hello, args=(i,q))
        processes.append(t)
        t.start()

    for one_process in processes:
        one_process.join()

    mylist = [ ]
    while not q.empty():
        mylist.append(q.get())

    print("Done!")
    print(len(mylist))
    print(mylist)

go()
