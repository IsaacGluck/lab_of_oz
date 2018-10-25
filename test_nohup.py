#!/usr/bin/env python
import time


time.sleep(60)
newFile = open('test_nohup.txt', 'w')
newFile.write('TESTING AFTER 1 MINUTE')
newFile.close()
