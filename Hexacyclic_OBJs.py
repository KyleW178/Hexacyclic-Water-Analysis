# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 15:51:35 2021

@author: Haworth_Lab2
"""
#A simplified version of water analysis, just for cyclic peptide structures

import os,math,csv,pytz,itertools
import numpy as np
import matplotlib.pyplot as plt
import datetime
import Hexacyclic_functions
#declear time
global starttime
starttime = str(datetime.datetime.now(pytz.timezone('America/Los_Angeles')))[0:16]


#============================
#Organize related data based on each structure
class structure():
    def __init__(self,folder:str = '',CA_list:dict = {},HB_pair:dict = {}):
        self.Master_folder = folder
        self.CA = CA_list
        self.HB_pair = HB_pair
#Zmat
    def atm_name_conv(self,origin_PDB:str = ''):
        #Open a pdb from PEPEDIT or other results and check the CA list/HB atm list is right
        
        pass


Main_folder = "C:\\People\\Kaichen\\cycGGGGGG\\"
AllG6 ='C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\ALL_G6'
AllA6 ='C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\ALL_A6'
AllV6 ='C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\V6'
All10 = "C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\ALL10"

debug = "C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\debug_test"
All_A6_1NMe = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\ALL_A6_1NMe'

ALL_A6_1_5NMe = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\ALL_A6_1_5NMe'
ALL_A6_1_2_5NMe = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\ALL_A6_1_2_5NMe'
ALL_A6_1_4_5NMe = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1,4,5NMe'
ALL_A6_1_2_4_5NMe = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\ALL_A6_1_2_4_5NMe'


ALL_A6_1_6NMe = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\ALL_A6_1_6NMe'
ALL_A6_1_3_6NMe = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\ALL_A6_1_3_6NMe'
ALL_A6_1_4_6NMe = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1,4,6NMe'
ALL_A6_1_3_4_6NMe = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\ALL_A6_1_3_4_6NMe'


#DAA
A6_1DAA = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1DAA'
ALL_A6_1NMe_1DAA = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1NMe_1DAA'
ALL_A6_1_5NMe_1DAA = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1,5NMe_1DAA'
ALL_A6_1_2_5NMe_1DAA = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1,2,5NMe_1DAA'
ALL_A6_1_4_5NMe_1DAA = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1,4,5NMe_1DAA'
ALL_A6_1_2_4_5NMe_1DAA = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1,2,4,5NMe_1DAA'


ALL_A6_1_6NMe_1DAA = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1_6NMe_1DAA'
ALL_A6_1_3_6NMe_1DAA ='C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1,3,6NMe_1DAA'
ALL_A6_1_4_6NMe_1DAA = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1,4,6NMe_1DAA'
ALL_A6_1_3_4_6NMe_1DAA = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-1,3,4,6NMe_1DAA'

#==================================
CA_G6 = {'backbone':[],'exposed':['C02','C04','C05','C07','C09','C11']}
CA_A6 = {'backbone':['C02','C05','C07','C10','C13','C16'],'exposed':['C03','C06','C09','C12','C15','C18']}
CA_V6 = {'backbone':['C02','C07','C11','C16','C21','C26'],'exposed':['C03','C04','C05','C08','C09','C10','C13','C14','C15','C18',
         'C19','C20','C23','C24','C25','C28','C29','C30']}

CA_A6_1NMe = {'backbone':['C02','C06','C08','C11','C14','C17'],'exposed':['C03','C04','C07','C10','C13','C16','C19']}

CA_A6_1_5NMe={'backbone':['C02','C06','C08','C11','C14','C18'],'exposed':['C03','C04','C07','C10','C13','C16','C17','C20']}
CA_A6_1_2_5NMe={'backbone':['C02','C06','C09','C12','C15','C19'],'exposed':['C03','C04','C07','C08','C11','C14','C17','C18','C21']}
CA_A6_1_4_5NMe={'backbone':['C02','C06','C08','C11','C15','C19'],'exposed':['C03','C04','C07','C10','C13','C14','C17','C18','C21']}

CA_A6_1_6NMe={'backbone':['C02','C06','C08','C11','C14','C17'],'exposed':['C03','C04','C07','C10','C13','C16','C19','C20']}
CA_A6_1_3_6NMe={'backbone':['C02','C06','C08','C12','C15','C18'],'exposed':['C03','C04','C07','C10','C11','C14','C17','C21','C20']}
CA_A6_1_4_6NMe ={'backbone':['C02','C06','C08','C11','C15','C18'],'exposed':['C03','C04','C07','C10','C13','C14','C17','C21','C20']}
CA_A6_1_3_4_6NMe={'backbone':['C02','C06','C08','C12','C16','C19'],'exposed':['C03','C04','C07','C10','C11','C14','C15','C18','C21','C22']}
CA_A6_1_2_4_5NMe ={'backbone':['C02','C06','C09','C12','C16','C20'],'exposed':['C03','C04','C07','C08','C11','C14','C15','C18','C19','C22']}

#DAA
CA_A6_1_6NMe_1DAA = {'backbone':['C02','C06','C08','C11','C14','C17'],'exposed':['C03','C04','C07','C10','C13','C16','C19','C20']}

#We can change this later.
HB_atm_pair_G6 = {'O01':'C01','O02':'C03','O03':'C06','O04':'C08','O05':'C10','O06':'C12',
               'N01':'H01','N02':'H02','N03':'H03','N04':'H04','N05':'H05','N06':'H06'}#FOR BASE TEMPLATE(all 6's) ONLY

HB_atm_pair_A6 = {'O01':'C01','O02':'C04','O03':'C08','O04':'C11','O05':'C14','O06':'C17',
               'N01':'H01','N02':'H02','N03':'H03','N04':'H04','N05':'H05','N06':'H06'}
HB_atm_pair_V6 = {'O01':'C01','O02':'C06','O03':'C12','O04':'C17','O05':'C22','O06':'C27',
               'N01':'H01','N02':'H02','N03':'H03','N04':'H04','N05':'H05','N06':'H06'}

HB_atm_pair_1NMe = {'O01':'C01','O02':'C05','O03':'C09','O04':'C12','O05':'C15','O06':'C18',
               'N02':'H01','N03':'H02','N04':'H03','N05':'H04','N06':'H05'}

HB_atm_pair_1_5NMe ={'O01':'C01','O02':'C05','O03':'C09','O04':'C12','O05':'C15','O06':'C19',
               'N02':'H01','N03':'H02','N04':'H03','N06':'H04'}
HB_atm_pair_1_2_5NMe ={'O01':'C01','O02':'C05','O03':'C10','O04':'C13','O05':'C16','O06':'C20',
                      'N03':'H01','N04':'H02','N06':'H03'}
HB_atm_pair_1_4_5NMe ={'O01':'C01','O02':'C05','O03':'C09','O04':'C12','O05':'C16','O06':'C20',
                      'N02':'H01','N03':'H02','N06':'H03'}

HB_atm_pair_1_6NMe ={'O01':'C01','O02':'C05','O03':'C09','O04':'C12','O05':'C15','O06':'C18',
               'N02':'H01','N03':'H02','N04':'H03','N05':'H04'}
HB_atm_pair_1_3_6NMe ={'O01':'C01','O02':'C05','O03':'C09','O04':'C13','O05':'C16','O06':'C19',
               'N02':'H01','N04':'H02','N05':'H03'}
HB_atm_pair_1_4_6NMe ={'O01':'C01','O02':'C05','O03':'C09','O04':'C12','O05':'C16','O06':'C19',
               'N02':'H01','N03':'H02','N05':'H03'}

HB_atm_pair_1_3_4_6NMe ={'O01':'C01','O02':'C05','O03':'C09','O04':'C13','O05':'C17','O06':'C20',
                         'N02':'H01','N05':'H02'}
HB_atm_pair_1_2_4_5NMe ={'O01':'C01','O02':'C05','O03':'C10','O04':'C13','O05':'C17','O06':'C21',
                         'N03':'H01','N06':'H02'}
#DAA
HB_atm_pair_1_6NMe_1DAA = {'O01':'C01','O02':'C05','O03':'C09','O04':'C12','O05':'C15','O06':'C18',
               'N02':'H01','N03':'H02','N04':'H03','N05':'H04'}
#
#HB_atm_pair_2NMe = {'O01':'C01','O02':'C04','O03':'C06','O04':'C08','O05':'C10','O06':'C12',
#               'N01':'H01','N03':'H02','N04':'H03','N05':'H04','N06':'H05'}


#=========================
G6 = structure(AllG6,CA_G6,HB_atm_pair_G6)
A6 = structure(AllA6,CA_A6,HB_atm_pair_A6)
V6 = structure(AllV6,CA_V6,HB_atm_pair_V6)
A6_1NMe = structure(All_A6_1NMe,CA_A6_1NMe,HB_atm_pair_1NMe)
A6_1_5NMe = structure(ALL_A6_1_5NMe,CA_A6_1_5NMe,HB_atm_pair_1_5NMe)
A6_1_2_5NMe = structure(ALL_A6_1_2_5NMe,CA_A6_1_2_5NMe,HB_atm_pair_1_2_5NMe)
A6_1_4_5NMe = structure(ALL_A6_1_4_5NMe,CA_A6_1_4_5NMe,HB_atm_pair_1_4_5NMe)
A6_1_2_4_5NMe = structure(ALL_A6_1_2_4_5NMe,CA_A6_1_2_4_5NMe,HB_atm_pair_1_2_4_5NMe)

A6_1_6NMe = structure(ALL_A6_1_6NMe,CA_A6_1_6NMe,HB_atm_pair_1_6NMe)
A6_1_3_6NMe = structure(ALL_A6_1_3_6NMe,CA_A6_1_3_6NMe,HB_atm_pair_1_3_6NMe)
A6_1_3_4_6NMe = structure(ALL_A6_1_3_4_6NMe,CA_A6_1_3_4_6NMe,HB_atm_pair_1_3_4_6NMe)
A6_1_4_6NMe = structure(ALL_A6_1_4_6NMe,CA_A6_1_4_6NMe,HB_atm_pair_1_4_6NMe)
#DAA corner
aA5 = structure(A6_1DAA,CA_A6,HB_atm_pair_A6)
A6_1NMe_1DAA = structure(ALL_A6_1NMe_1DAA,CA_A6_1NMe,HB_atm_pair_1NMe)
A6_1_5NMe_1DAA = structure(ALL_A6_1_5NMe_1DAA,CA_A6_1_5NMe,HB_atm_pair_1_5NMe)
A6_1_2_5NMe_1DAA = structure(ALL_A6_1_2_5NMe_1DAA,CA_A6_1_2_5NMe,HB_atm_pair_1_2_5NMe)
A6_1_4_5NMe_1DAA = structure(ALL_A6_1_4_5NMe_1DAA,CA_A6_1_4_5NMe,HB_atm_pair_1_4_5NMe)
A6_1_2_4_5NMe_1DAA = structure(ALL_A6_1_2_4_5NMe_1DAA,CA_A6_1_2_4_5NMe,HB_atm_pair_1_2_4_5NMe)

A6_1_6NMe_1DAA = structure(ALL_A6_1_6NMe_1DAA,CA_A6_1_6NMe_1DAA,HB_atm_pair_1_6NMe_1DAA)
A6_1_3_6NMe_1DAA = structure(ALL_A6_1_3_6NMe_1DAA,CA_A6_1_3_6NMe,HB_atm_pair_1_3_6NMe)
A6_1_4_6NMe_1DAA =structure(ALL_A6_1_4_6NMe_1DAA,CA_A6_1_4_6NMe,HB_atm_pair_1_4_6NMe)
A6_1_3_4_6NMe_1DAA = structure(ALL_A6_1_3_4_6NMe_1DAA,CA_A6_1_3_4_6NMe,HB_atm_pair_1_3_4_6NMe)

#1DAA-NMe rotation
#A6_2NMe_2DAA = structure(ALL_A6_2NMe_2DAA,CA_A6_2NMe,HB_atm_pair_2NMe)


#===========================






                    
if __name__ == '__main__':
#        
#    #print(get_report("C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\debug_test\\Pose_3381"))
#    
#    
#    output_report(A6_1_5NMe_1DAA,opt = 'csv',graph = 'N')
    #print(get_report(A6,debug+'\\Pose_1010'))


    new_lims = permutation([4,8,4,0.05],{1:[6,8],2:[3,4]})

    graph_output(A6,new_lims)
    
