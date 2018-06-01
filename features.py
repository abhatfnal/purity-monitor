
def ProgressBar(it, counts):
    width = counts
    if(it==0):
        sys.stdout.write("[%s]" % (" " * width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width+1))
        sys.stdout.write("-")
    else:
        sys.stdout.write("-")
        sys.stdout.flush()
    if(it==counts):
        sys.stdout.write("\n")
