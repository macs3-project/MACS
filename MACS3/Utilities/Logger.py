# Logger.py to set time and memory monitoring to logging
import logging
import resource
import time
import os
import sys

class MemoryLogger(logging.Logger):
    def __init__(self, name, level=logging.NOTSET):
        super().__init__(name, level)
    
    def _log(self, level, msg, args, exc_info=None, extra=None, stack_info=False):
        mem_usage = self.get_memory_usage()
        super()._log(level, f"[{mem_usage} MB] {msg}", args, exc_info, extra, stack_info)

    @staticmethod
    def get_memory_usage():
        mem_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if os.name == 'posix' and os.uname().sysname == 'Darwin':
            # macOS
            mem_usage = mem_usage / 1024  # Convert to kilobytes
        return int( mem_usage / 1024 ) # Convert to MB

logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

logging.setLoggerClass(MemoryLogger)