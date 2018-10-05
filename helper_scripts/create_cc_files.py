
import sys, os
import shutil, glob

# given a target directory: tfolder
# containing some .hpp files, this scripts creates
# a tfolder/src directory where it puts all .cc files
# and then puts and include inside these .cc files to include .hpp ones

def createSrc(dest):
    # check if src exists
    srcd = os.getcwd()+"/src"
    bo = os.path.isdir(srcd)
    # if exists, wipe it
    if bo:
        shutil.rmtree(srcd)

    mylist = glob.glob("*.hpp")
    print mylist
    os.mkdir(dest+"/src")
    for ifile in mylist:
        hppname = ifile.replace('.hpp','')
        ccname = hppname + '.cc'
        os.system('touch ./src/' + ccname)
        incstr = '#include "../'+hppname+'.hpp"'
        out = open("./src/"+ccname, "w")
        out.write(incstr)
        out.close()
  
def createCCFiles(path, cwd):
    os.chdir(path)
    thisF = os.getcwd()
    # get list of dirs here
    dirs_here = filter(os.path.isdir, os.listdir(os.curdir))
    # remove the 'src' from the list if there is one
    dirs2 = [i for i in dirs_here if i!="src" and i!="tinympl"]

    # create for current src
    root = os.getcwd()
    createSrc(root)
    # loop over subfolders
    if dirs2:
        for jj in dirs2:
            createCCFiles(jj,thisF)
    else: # otherwise, if no subfolders, create src and put .cc files
        root = os.getcwd()
        print ("we are in: " + root)
        createSrc(root)
    os.chdir(cwd)
    
# Check current working directory.
cwd = os.getcwd()
print "Current working directory %s" % cwd

tFolder = cwd+"/../packages"
pn = ['core', 'solvers', 'svd', 'ode', 'rom']
for i in pn:
    path = tFolder+"/"+i+"/src"
    print ("\nProcessing package **" + i + "** in " + path)
    thisf = os.getcwd()
    # creat 
    createCCFiles(path, thisf)
