import unittest
import os
import tempfile
import shutil
from pkg_resources import resource_filename
import sys
import subprocess
import glob
import scipy as sp
if (sys.version_info > (3, 0)):
    # Python 3 code in this block
    import configparser as ConfigParser
else:
    import ConfigParser

class TestCor(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._branchFiles = tempfile.mkdtemp()+"/"

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls._branchFiles):
            shutil.rmtree(cls._branchFiles, ignore_errors=True)

    def test_cor(self):

        self._test = True
        self.send_requirements()
        self.send_submit_mocks()

        if self._test:
            self.remove_folder()

        return
    def produce_folder(self):
        """
            Create the necessary folders
        """

        print("\n")
        lst_fold = ['/Products/']

        for fold in lst_fold:
            if not os.path.isdir(self._branchFiles+fold):
                os.mkdir(self._branchFiles+fold)

        return
    def remove_folder(self):
        """
            Remove the produced folders
        """

        print("\n")
        shutil.rmtree(self._branchFiles, ignore_errors=True)

        return
    def load_requirements(self):

        req = {}
        os.environ['SACLAYMOCKS_BASE'] = resource_filename('SaclayMocks', '/../../')
        path = os.path.expandvars('$SACLAYMOCKS_BASE/requirements.txt')
        with open(path,'r') as f:
            for l in f:
                l = l.replace('\n','').replace('==',' ').replace('>=',' ').split()
                if 'git+https://github.com/' in l[0]:
                    l[0] = l[0].split('.git')[-2].split('/')[-1]
                if len(l)==1:
                    req[l[0]] = 0
                elif len(l)==2:
                    req[l[0]] = l[1]
                else:
                    self.assertTrue(False,"requirements.txt attribute is not valid: {}".format(str(l)))

        return req



    def send_requirements(self):

        print("\n")
        req = self.load_requirements()
        for req_lib, req_ver in req.items():
            try:
                local_ver = __import__(req_lib).__version__
                if local_ver!=req_ver:
                    print("WARNING: The local version of {}: {} is different from the required version: {}".format(req_lib,local_ver,req_ver))
            except ModuleNotFoundError:
                print("WARNING: Module {} has some errors".format(req_lib))
            except ImportError:
                print("WARNING: Module {} can't be found".format(req_lib))
            except AttributeError:
                print("WARNING: Module {} has no version".format(req_lib))

        return
    def send_submit_mocks(self):

        ###
        print("\n")
        cmd = 'submit_mocks.py'
        cmd += ' --mock-dir '+self._branchFiles+'/Products/'
        cmd += ' --cori-nodes False'
        cmd += ' --box-size 128'
        cmd += ' --chunk-id 1'
        cmd += ' --seed 42'
        print(cmd)
        subprocess.call(cmd, shell=True)

        ###
        print("\n")
        cmd = self._branchFiles+'/Products/mock_0/output/runs/submit.sh'
        print(cmd)
        subprocess.call(cmd, shell=True)

        ###
        print("\n")
        tl = ''
        for i in range(10):
            tl += '/*/'
            fs = glob.glob(self._branchFiles+'/Products/{}/*.log*'.format(tl))
            if len(fs)==0: continue
            print(fs)
            print("\n",i,"\n")
            for f in fs:
                with open(f) as tf:
                    print(tf.read())

        ###
        print("\n")
        tl = ''
        for i in range(10):
            tl += '/*/'
            fs = glob.glob(self._branchFiles+'/Products/{}/*.log*'.format(tl))
            if len(fs)==0: continue
            print(fs)
            print("\n",i,"\n")
            for f in fs:
                with open(f) as tf:
                    tff = sp.hstack([ l.split() for l in tf ])
                    self.assertFalse('Abort!' in tff,"ERROR: 'Abort!' in {}".format(f))
                    self.assertFalse('Error' in tff,"ERROR: 'Error' in {}".format(f))

        ### Test
        #if self._test:
        #    path1 = self._masterFiles + '/Products/'
        #    path2 = self._branchFiles + '/Products/'
        #    self.compare_fits(path1,path2,"submit_mocks.py")

        return

if __name__ == '__main__':
    unittest.main()
