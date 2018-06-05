import sys

def ProgressBar(it, counts):
    width = counts
    if(it==1):
        sys.stdout.write(" | [%s]" % (" " * width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (width+1))
        sys.stdout.write("-")
    else:
        sys.stdout.write("-")
        sys.stdout.flush()
    if(it==counts):
        sys.stdout.write("\n")

def ElapsedTime():
    print " | Time elapsed:          ", time.clock() , "sec"
