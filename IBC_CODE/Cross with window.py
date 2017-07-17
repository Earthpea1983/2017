# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 14:00:08 2017

@author: k-hlm
"""

import os
import re

import pandas as pd
import win32ui


class Cross_cargo_list:
    def __init__(self):
        self.openfile_IBC()
        self.openfile_cargolist()
        self.cross_list=pd.DataFrame([],columns=self.ibc.columns)
        self.compare()
        self.save2excel()
        self.Texproof()
        self.Aexproof()
    
    def openfile_IBC(self): #open IBC
        dlg=win32ui.CreateFileDialog(1)
        dlg.SetOFNInitialDir(os.getcwd())
        dlg.DoModal()
        filename=dlg.GetPathName()
        try:
            self.ibc=pd.read_excel(filename)
        except:
            print("无法打开文件<IBC_cargo_list.xlsx>!")
    
    def openfile_cargolist(self):
        dlg=win32ui.CreateFileDialog(1)
        dlg.SetOFNInitialDir(os.getcwd())
        dlg.DoModal()
        filename=dlg.GetPathName()        
        try:
            self.cargo_list=pd.read_excel(filename)
        except:
            print("无法打开文件cargo_list")            
            
    def compare(self): #cross cargo list and IBC
        for j in self.cargo_list.Name:
            pattern=r"{}".format(j)
            for i in range(len(self.ibc)):
                if re.search(pattern,self.ibc.ix[i,0]):
                    self.cross_list=self.cross_list.append(self.ibc.ix[i,:])
        
    def save2excel(self): #save file to excel
        self.cross_list.to_excel("Cross_result.xlsx")
        print("结果存储在<Cross_result.xlsx>!")
        
    def Texproof(self): #find the max temperature class of ex
        for j in range(7,0,-1):
            pattern=r"{}".format(j)
            indT=self.cross_list.index
            for i in indT:
                if str(self.cross_list.ix[i,'Temperature_Classes'])=='nan':
                    continue
                if re.search(pattern,str(self.cross_list.ix[i,'Temperature_Classes'])):
                    return print("最大温度等级{}".format("T"+str(j))+": {}".format(self.cross_list.ix[i,'Product_Name']))
        print("货品无温度等级要求！")
        
    def Aexproof(self): # fine the max apparatus class of ex
        Apparatus=['C','B','A']
        for j in Apparatus:
            pattern=r"{}".format(j)
            indA=self.cross_list.index
            for i in indA:
                if str(self.cross_list.ix[i,'Apparatus_group'])=='nan':
                    continue
                if re.search(pattern,str(self.cross_list.ix[i,'Apparatus_group'])):
                    return print("最大间隙等级{}".format("II"+str(j))+": {}".format(self.cross_list.ix[i,'Product_Name']))                        
        print("货品无间隙等级要求！")
        
        
my=Cross_cargo_list()                    
        