#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:09:53 2021

@author: sohailnizam
"""

#-------------------------------------------------#
# Scrape cvR2 and CI data from slapnap html files #
#-------------------------------------------------#


#import necessary libraries
import os
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
from glob import glob

#the bnabs which have directories in ic50ic80/docker_output 
#bnabs = ["vrc01", "10-1074", "10-996", "2f5", "2g12", "35o22",
#         "3bnc117", "4e10", "8anc195", "b12", "ch01", "hj16",
#         "nih45-46", "pg16", "pg9", "pgdm1400", "pgt121", 
#         "pgt128", "pgt135", "pgt145", "pgt151", "vrc-ch31",
#        "vrc-pg04", "vrc01", "vrc03", "vrc07", "vrc26.08",
#        "vrc26.25", "vrc29.03", "vrc34.01"]

bnabs = os.listdir("./docker_output")


def scrape_tables(bnabs):
    
    '''
    This function iterates over the directories of the bnabs,
    extracts the cross validated R2 as well as the CI lb and ubs
    for both ic50 and ic80. These values are stored in a pandas df.
    
    Input: list of bnab name strings
    Output: pandas df with 2 rows per bnab (1 for ic50, 1 for ic80)
            and 4 columns (bnab, method, r2, cil, ciu)
    '''
    #init empty list to hold scraped r2s, lbs, and ubs
    cvR2_data = []

    
    for bnab in bnabs:
        
        #assume ic50ic80 is the current working dir
        path = "./docker_output/"
        path += bnab
        #report_path = path + "/report_" + bnab.upper() + "_" + date + ".html"
        
        
        #if the container failed, the dir does not have a report
        try:
            #pick out the file that starts with 'report' ends with '.html'
            #this method is not date specific
            report_path = glob(path + "/report*.html")[0]
            report = open(report_path).read()
            
        except IndexError:
            
            break
        
        #scrape the html with beautiful soup
        soup = BeautifulSoup(report)
        #extract the second table in the html
        cvR2_table = soup.findAll("table")[1]
        
        #init a list to hold the data for this bnab
        output_rows = []
        #two rows not including header: ic50 and ic80
        for table_row in cvR2_table.findAll("tr")[1:]:
            
            #three columns not including label
            columns = table_row.findAll('td')
            output_row = []
            for column in columns[1:]:
            
                output_row.append(column.text)
                
            output_rows.append(output_row)
        
        #add the rows for this bnab to the running list
        cvR2_data += output_rows
        #print(bnab + " done")
    
    
    #cast the data list of lists to pandas df
    cvR2_dataframe = pd.DataFrame(cvR2_data, columns = ["r2", "cil", "ciu"])
    
    
    #add nab and ic columns
    cvR2_dataframe.insert(loc=0, column='bnab', value=np.repeat(bnabs, 2))
    cvR2_dataframe.insert(loc=1, column='method', value=['ic50', 'ic80']*len(bnabs))


    #return sorted to ic50 and ic80 are separate
    return(cvR2_dataframe.sort_values(by = ['method', 'bnab']))
        

#call the function on the bnabs
result_df = scrape_tables(bnabs)
        
#write to csv
result_df.to_csv("./Python/cvr2_df.csv")
 

    
    
    