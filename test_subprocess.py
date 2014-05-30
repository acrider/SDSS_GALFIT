def run_galfit():

    p = subprocess.Popen(['/Users/acrider/galfit/galfit', 'galfit.feedme'], cwd='/Users/acrider/Desktop/my-example',
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        print line,
    retval = p.wait()
    
    return 

def run_ls():

    p = subprocess.Popen('ls', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        print line,
    retval = p.wait()
    
    return
    
#----

import subprocess

run_ls()

run_galfit()
