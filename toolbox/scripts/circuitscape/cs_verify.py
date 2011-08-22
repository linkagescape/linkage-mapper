##
## Circuitscape (C) 2008, 2009, 2010, Brad McRae and Viral B. Shah. 
##
## $Id: cs_verify.py 695 2009-11-04 18:43:52Z mcrae $
##

import imp, os, sys

import unittest

from cs_util import *
from cs_compute import *
import numpy as np
np.seterr(all = 'ignore')

def approxEqual(a, b):

    m = a.shape[0]
    n = a.shape[1]

    for i in range(0, m):
        for j in range(0, n):
            if (a[i,j] != b[i,j]):
                if (abs(a[i,j] - b[i,j]) > 1e-6):
                    return False
    return True

def cs_verifyall():
    suite = unittest.TestLoader().loadTestsFromTestCase(cs_verify)
    testResult = unittest.TextTestRunner(verbosity=0).run(suite)
    return testResult
    unittest.main()

def test_sg(ut, test_name):
    #print test_name
    configFile='.//verify//config_files//' + test_name + '.ini'
    cs = cs_compute(configFile, None)
    resistances_computed,solver_failed = cs.compute()
    resistances_saved=loadtxt('.//verify//baseline_results//' + test_name + '_resistances.txt') 
    
    ut.assertEquals (approxEqual(resistances_saved, resistances_computed), True)    


    if test_name!='sgVerify12': #This module tests the shortcut which is only calculated when maps are not written.
        current_map_1_2_computed=reader('.//verify//output//' + test_name + '_curmap_1_2.asc', 'float64') 
        current_map_1_2_saved=reader('.//verify//baseline_results//' + test_name + '_curmap_1_2.asc', 'float64') 
        cum_current_map_computed=reader('.//verify//output//' + test_name + '_cum_curmap.asc', 'float64') 
        cum_current_map_saved=reader('.//verify//baseline_results//' + test_name + '_cum_curmap.asc', 'float64')     
        voltage_map_1_2_computed=reader('.//verify//output//' + test_name + '_voltmap_1_2.asc', 'float64') 
        voltage_map_1_2_saved=reader('.//verify//baseline_results//' + test_name + '_voltmap_1_2.asc', 'float64') 
            
        ut.assertEquals (approxEqual(current_map_1_2_saved, current_map_1_2_computed), True)
        ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
        ut.assertEquals (approxEqual(voltage_map_1_2_saved, voltage_map_1_2_computed), True)
    

def test_one_to_all(ut, test_name):
    #print test_name
    configFile='.//verify//config_files//' + test_name + '.ini'
    cs = cs_compute(configFile, None)
    resistances_computed,solver_failed = cs.compute()

    resistances_saved=loadtxt('.//verify//baseline_results//' + test_name + '_resistances.txt') 

    current_map_1_computed=reader('.//verify//output//' + test_name + '_curmap_1.asc', 'float64') 
    current_map_1_saved=reader('.//verify//baseline_results//' + test_name + '_curmap_1.asc', 'float64') 

    cum_current_map_computed=reader('.//verify//output//' + test_name + '_cum_curmap.asc', 'float64') 
    cum_current_map_saved=reader('.//verify//baseline_results//' + test_name + '_cum_curmap.asc', 'float64') 

    voltage_map_1_computed=reader('.//verify//output//' + test_name + '_voltmap_1.asc', 'float64') 
    voltage_map_1_saved=reader('.//verify//baseline_results//' + test_name + '_voltmap_1.asc', 'float64') 
    
    ut.assertEquals (approxEqual(resistances_saved, resistances_computed), True)
    ut.assertEquals (approxEqual(current_map_1_saved, current_map_1_computed), True)
    ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
    ut.assertEquals (approxEqual(voltage_map_1_saved, voltage_map_1_computed), True)

def test_all_to_one(ut, test_name):
    #print test_name
    configFile='.//verify//config_files//' + test_name + '.ini'
    cs = cs_compute(configFile, None)
    resistances_computed,solver_failed = cs.compute()

    current_map_1_computed=reader('.//verify//output//' + test_name + '_curmap_1.asc', 'float64') 
    current_map_1_saved=reader('.//verify//baseline_results//' + test_name + '_curmap_1.asc', 'float64') 

    cum_current_map_computed=reader('.//verify//output//' + test_name + '_cum_curmap.asc', 'float64') 
    cum_current_map_saved=reader('.//verify//baseline_results//' + test_name + '_cum_curmap.asc', 'float64') 

    voltage_map_1_computed=reader('.//verify//output//' + test_name + '_voltmap_1.asc', 'float64') 
    voltage_map_1_saved=reader('.//verify//baseline_results//' + test_name + '_voltmap_1.asc', 'float64') 
    
    ut.assertEquals (approxEqual(current_map_1_saved, current_map_1_computed), True)
    ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
    ut.assertEquals (approxEqual(voltage_map_1_saved, voltage_map_1_computed), True)


        
def test_mg(ut, test_name):
    #print test_name
    configFile='.//verify//config_files//' + test_name + '.ini'
    cs = cs_compute(configFile, None)
    voltages = cs.compute()
   
    cum_current_map_computed=reader('.//verify//output//' + test_name + '_curmap.asc', 'float64') 
    cum_current_map_saved=reader('.//verify//baseline_results//' + test_name + '_curmap.asc', 'float64') 

    voltage_map_computed=reader('.//verify//output//' + test_name + '_voltmap.asc', 'float64') 
    voltage_map_saved=reader('.//verify//baseline_results//' + test_name + '_voltmap.asc', 'float64') 
    
    ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
    ut.assertEquals (approxEqual(voltage_map_saved, voltage_map_computed), True)
        
class cs_verify(unittest.TestCase):
    def test_single_ground_all_pairs_resistances_1(self):
        test_sg(self, 'sgVerify1') 

    def test_single_ground_all_pairs_resistances_2(self):
        test_sg(self, 'sgVerify2') 

    def test_single_ground_all_pairs_resistances_3(self):
        test_sg(self, 'sgVerify3') 

    def test_single_ground_all_pairs_resistances_4(self):
        test_sg(self, 'sgVerify4') 

    def test_single_ground_all_pairs_resistances_5(self):
        test_sg(self, 'sgVerify5') 

    def test_single_ground_all_pairs_resistances_6(self):
        test_sg(self, 'sgVerify6') 

    def test_single_ground_all_pairs_resistances_7(self):
        test_sg(self, 'sgVerify7') 

    def test_single_ground_all_pairs_resistances_8(self):
        test_sg(self, 'sgVerify8') 

    def test_single_ground_all_pairs_resistances_9(self):
        test_sg(self, 'sgVerify9') 

    def test_single_ground_all_pairs_resistances_10(self):
        test_sg(self, 'sgVerify10')

    def test_single_ground_all_pairs_resistances_11(self):
        test_sg(self, 'sgVerify11') 
    
    def test_single_ground_all_pairs_resistances_12(self):
        test_sg(self, 'sgVerify12') 
        
    def test_single_ground_all_pairs_resistances_13(self):
        test_sg(self, 'sgVerify13')         
        
    def test_multiple_ground_module_1(self):
        test_mg(self, 'mgVerify1') 
      
    def test_multiple_ground_module_2(self):
        test_mg(self, 'mgVerify2') 

    def test_multiple_ground_module_3(self):
        test_mg(self, 'mgVerify3') 
        
    def test_multiple_ground_module_4(self):
        test_mg(self, 'mgVerify4')         

    def test_multiple_ground_module_5(self):
        test_mg(self, 'mgVerify5')         
  
    def test_one_to_all_module_1(self):
        test_one_to_all(self, 'oneToAllVerify1') 
  
    def test_one_to_all_module_2(self):
        test_one_to_all(self, 'oneToAllVerify2') 

    def test_one_to_all_module_3(self):
        test_one_to_all(self, 'oneToAllVerify3') 

    def test_one_to_all_module_4(self):
        test_one_to_all(self, 'oneToAllVerify4') 

    def test_one_to_all_module_5(self):
        test_one_to_all(self, 'oneToAllVerify5') 
  
    def test_one_to_all_module_6(self):
        test_one_to_all(self, 'oneToAllVerify6') 

    def test_one_to_all_module_7(self):
        test_one_to_all(self, 'oneToAllVerify7') 

    def test_one_to_all_module_8(self):
        test_one_to_all(self, 'oneToAllVerify8') 

    def test_one_to_all_module_9(self):
        test_one_to_all(self, 'oneToAllVerify9') 
  
    def test_one_to_all_module_10(self):
        test_one_to_all(self, 'oneToAllVerify10') 

    def test_one_to_all_module_11(self):
        test_one_to_all(self, 'oneToAllVerify11') 

    def test_one_to_all_module_12(self):
        test_one_to_all(self, 'oneToAllVerify12') 
        
    def test_all_to_one_module_1(self):
        test_all_to_one(self, 'allToOneVerify1') 
  
    def test_all_to_one_module_2(self):
        test_all_to_one(self, 'allToOneVerify2') 

    def test_all_to_one_module_3(self):
        test_all_to_one(self, 'allToOneVerify3') 

    def test_all_to_one_module_4(self):
        test_all_to_one(self, 'allToOneVerify4') 

    def test_all_to_one_module_5(self):
        test_all_to_one(self, 'allToOneVerify5') 
  
    def test_all_to_one_module_6(self):
        test_all_to_one(self, 'allToOneVerify6') 

    def test_all_to_one_module_7(self):
        test_all_to_one(self, 'allToOneVerify7') 

    def test_all_to_one_module_8(self):
        test_all_to_one(self, 'allToOneVerify8') 

    def test_all_to_one_module_9(self):
        test_all_to_one(self, 'allToOneVerify9') 
  
    def test_all_to_one_module_10(self):
        test_all_to_one(self, 'allToOneVerify10') 

    def test_all_to_one_module_11(self):
        test_all_to_one(self, 'allToOneVerify11')   
        
    def test_all_to_one_module_12(self):
        test_all_to_one(self, 'allToOneVerify12')           
if __name__ == '__main__':
    unittest.main()
